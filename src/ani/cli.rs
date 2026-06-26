//! Command-line argument parsing and the `--help` text.

use std::{
    env, fs, io, path,
    path::{Path, PathBuf},
    process,
};

use crate::ani::{
    default_max_shard_minimizers_for_runtime, describe_field_parsing_error, is_stdin_path,
    load_params_file, validate_fragment_length, validate_kmer_size, validate_mash_confidence,
    validate_max_shard_minimizers, validate_min_identity, validate_window_size, FastaInput,
    IndexBuildMode, ParamsFileConfig, DEFAULT_FRAGMENT_LENGTH, DEFAULT_FRAGMENT_STRIDE,
    DEFAULT_FREQ_THRESHOLD_PERCENT, DEFAULT_KMER_SIZE, DEFAULT_MASH_CONFIDENCE,
    DEFAULT_MAX_SHARD_MINIMIZERS, DEFAULT_MIN_FRAGMENT_LENGTH, DEFAULT_MIN_PERCENT_IDENTITY,
    DEFAULT_SPLIT_N_RUN, DEFAULT_WINDOW_SIZE,
};

/// Parsed command-line arguments.
pub(crate) struct CliArgs {
    pub(crate) references: Vec<FastaInput>,
    pub(crate) queries: Vec<FastaInput>,
    pub(crate) sketch_path: Option<PathBuf>,
    pub(crate) tmp_dir: Option<PathBuf>,
    pub(crate) out_path: Option<PathBuf>,
    pub(crate) mapping_stats_path: Option<PathBuf>,
    pub(crate) bgzip: bool,
    pub(crate) emit_header: bool,
    pub(crate) verbose: bool,
    pub(crate) threads: usize,
    pub(crate) freq_threshold_percent: f64,
    pub(crate) minmer_count: Option<usize>,
    pub(crate) kmer_size: usize,
    pub(crate) window_size: usize,
    pub(crate) fragment_length: u32,
    pub(crate) fragment_stride: u32,
    pub(crate) min_fragment_length: u32,
    pub(crate) min_identity: f64,
    pub(crate) mash_confidence: f64,
    /// Minimum run-length of ambiguous `N` bases that splits a contig; `0` disables.
    pub(crate) split_n_run: usize,
    pub(crate) max_memory_bytes: Option<u64>,
    pub(crate) max_shard_minimizers: usize,
    pub(crate) max_concurrent_shards: Option<usize>,
    pub(crate) index_build_mode: IndexBuildMode,
}

fn usage() -> &'static str {
    "usage: fasterANI (--reference <ref.fa> | --reference-list <refs.txt>)... \
[(--query <query.fa> | --query-list <queries.txt>)...] [options]

Inputs:
  --params-file <path>         Load runtime parameters from a TOML file.
                                Values in the file override defaults; CLI arguments override
                                  file values. Relative paths in the file are resolved relative
                                  to the TOML file's directory.
  --reference <path>            Reference FASTA (optionally gzip-compressed).
  --reference-list <path>       File of reference FASTA paths, one per line.
  --reference-sketch <prefix>   Build/reuse an on-disk reference sketch at this prefix.
                                Query inputs optional (build-only when omitted).
  --query <path>                Query FASTA (repeatable). Use `-` to read one query from
                                  stdin (optionally gzip-compressed).
  --query-list <path>           File of query FASTA paths, one per line.
  --query-name <label>          Display/path label for a stdin query (with `--query -`).

Output:
  --out <path>                  Write results TSV here (default: stdout).
  --header                      Prepend a column-name header row to the results TSV
                                  (default: off, for FastANI/script compatibility).
  --mapping-stats <path>        Write a per-fragment mapping-stats TSV (always headered).
  --verbose                     Print PROGRESS/diagnostics to stderr (default: off).

  Results columns (tab-separated):
    query_file  reference_file  ani  shared_fragment_equivalents  total_fragment_equivalents
  where ANI is a percent, and the last two are fractional fragment counts
  (aligned bases / fragment-length), so they may be non-integer.

Seeding (minimizer sketch; applies to both references and queries):
  --kmer-size <n>               K-mer size for minimizers (default 16).
  --window-size <n>             Minimizer window size (default 24).
  --minmer-count <n>            Keep only the n smallest-hash minimizers ('minmers') per query
                                  fragment for candidate scoring (default behavior uses all).
  --freq-threshold-percent <p>  Ignore reference minimizers occurring in more than p%
                                  of reference positions; 0..100, 0 disables (default 0).

Fragmenting (how each query contig is cut into fragments):
  --fragment-length <bp>        Query fragment length (default 3000).
  --fragment-stride <bp>        Step between fragment starts; <= fragment-length
                                  (default: equal to fragment-length, i.e. non-overlapping).
  --min-fragment-length <bp>    Keep trailing fragments at least this long
                                  (alias: --min-fraglen; default: fragment-length).
  --split-N <bp>                Split contigs at runs of >= this many ambiguous (N)
                                  bases; 0 disables splitting (alias: --split-n; default 0).

Fragment mapping (thresholds applied to each individual fragment alignment):
  --min-identity <percent>      Minimum identity of a single query-fragment-to-reference
                                  alignment for that fragment to count toward ANI;
                                  0..100 (default 80). (Not a threshold on the final ANI.)
  --mash-confidence <fraction>  Confidence level for the Mash-distance prefilter that
                                  selects which reference regions each fragment is scored
                                  against; higher = stricter; 0..1 (default 0.9).

  Per genome pair, fasterANI keeps only reciprocal-best fragment mappings and reports
  ANI as the length-weighted mean of those retained fragments' identities.

Sketch database / sharding:
  --bgzip                       Enable if reference sketch input is bgzip-compressed.
  --max-shard-minimizers <n>    Maximum estimated reference minimizers per shard
                                  (default is 500_000_000, producing ~10 GiB shards).
  --max-concurrent-shards <n>    Maximum number of shards to query concurrently
  --index-build-mode <mode>     auto | hash | partitioned (default auto).

Resources:
  --threads <n>                 Worker threads, >= 1 (default 1).
  --max-memory-gb <gb>          Soft memory ceiling in GB (default: unlimited).
  --tmp <dir>                   Directory for temporary shard files.
  -h, --help                    Show help.
  -v, --version                 Show version."
}

fn validate_and_read_path_list_from_base(
    list_path: &str,
    base_dir: Option<&Path>,
) -> io::Result<(PathBuf, Vec<String>)> {
    let absolute_path = resolve_path_from_base(list_path, base_dir)?;
    let contents = fs::read_to_string(&absolute_path)?;
    let mut valid_paths = Vec::new();
    let list_base_dir = if base_dir.is_some() {
        absolute_path.parent()
    } else {
        None
    };

    for line in contents.lines() {
        let path = line.trim();
        if path.is_empty() || path.starts_with('#') {
            continue;
        }

        let absolute_item_path = resolve_path_from_base(path, list_base_dir)?;

        match fs::metadata(&absolute_item_path) {
            Ok(meta) if meta.is_file() => {
                valid_paths.push(absolute_item_path.to_string_lossy().into_owned());
            }
            Ok(_) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("Path in list is not a file: {}", path),
                ));
            }
            Err(e) => {
                return Err(io::Error::new(
                    e.kind(),
                    format!(
                        "Cannot access path '{}' from list: {}\n{}",
                        path, list_path, e
                    ),
                ));
            }
        }
    }
    Ok((absolute_path, valid_paths))
}

fn resolve_path_from_base(value: &str, base_dir: Option<&Path>) -> io::Result<PathBuf> {
    let path = PathBuf::from(value);
    if path.is_absolute() {
        return path::absolute(path);
    }
    match base_dir {
        Some(base_dir) => path::absolute(base_dir.join(path)),
        None => path::absolute(path),
    }
}

#[derive(Default)]
struct RuntimeStartupOutput {
    params_file: Option<StartupValue>,
    reference_files: Vec<StartupValue>,
    reference_lists: Vec<StartupValue>,
    query_files: Vec<StartupValue>,
    query_lists: Vec<StartupValue>,
    query_name: Option<StartupValue>,
}

#[derive(Clone, Copy, Eq, PartialEq)]
enum ParameterSource {
    ParamsFile,
    Cli,
}

impl ParameterSource {
    fn label(self) -> &'static str {
        match self {
            Self::ParamsFile => "params file",
            Self::Cli => "CLI",
        }
    }
}

struct StartupValue {
    value: String,
    source: ParameterSource,
}

impl StartupValue {
    fn new(value: impl Into<String>, source: ParameterSource) -> Self {
        Self {
            value: value.into(),
            source,
        }
    }
}

impl RuntimeStartupOutput {
    fn emit(&self, args: &CliArgs, sources: &ParameterSources) {
        let mut entries: Vec<String> = Vec::new();

        push_toml_string(&mut entries, "params_file", self.params_file.as_ref());
        push_toml_array(&mut entries, "reference_files", &self.reference_files);
        push_toml_array(&mut entries, "reference_lists", &self.reference_lists);
        push_toml_array(&mut entries, "query_files", &self.query_files);
        push_toml_array(&mut entries, "query_lists", &self.query_lists);
        push_toml_string(&mut entries, "query_name", self.query_name.as_ref());
        push_toml_path(
            &mut entries,
            "reference_sketch",
            args.sketch_path.as_ref(),
            sources.reference_sketch,
        );
        push_toml_path(&mut entries, "tmp", args.tmp_dir.as_ref(), sources.tmp);
        push_toml_path(&mut entries, "out", args.out_path.as_ref(), sources.out);
        push_toml_path(
            &mut entries,
            "mapping_stats",
            args.mapping_stats_path.as_ref(),
            sources.mapping_stats,
        );
        push_toml_bool(&mut entries, "bgzip", args.bgzip, sources.bgzip);
        push_toml_bool(&mut entries, "header", args.emit_header, sources.header);
        push_toml_bool(&mut entries, "verbose", args.verbose, sources.verbose);

        if sources.threads.is_some() && args.threads != 1 {
            push_toml_number(&mut entries, "threads", args.threads, sources.threads);
        }
        if sources.freq_threshold_percent.is_some()
            && args.freq_threshold_percent != DEFAULT_FREQ_THRESHOLD_PERCENT
        {
            push_toml_number(
                &mut entries,
                "freq_threshold_percent",
                args.freq_threshold_percent,
                sources.freq_threshold_percent,
            );
        }
        if let Some(minmer_count) = args.minmer_count {
            push_toml_number(
                &mut entries,
                "minmer_count",
                minmer_count,
                sources.minmer_count,
            );
        }
        if sources.kmer_size.is_some() && args.kmer_size != DEFAULT_KMER_SIZE {
            push_toml_number(&mut entries, "kmer_size", args.kmer_size, sources.kmer_size);
        }
        if sources.window_size.is_some() && args.window_size != DEFAULT_WINDOW_SIZE {
            push_toml_number(
                &mut entries,
                "window_size",
                args.window_size,
                sources.window_size,
            );
        }
        if sources.fragment_length.is_some() && args.fragment_length != DEFAULT_FRAGMENT_LENGTH {
            push_toml_number(
                &mut entries,
                "fragment_length",
                args.fragment_length,
                sources.fragment_length,
            );
        }
        if sources.fragment_stride.is_some() && args.fragment_stride != args.fragment_length {
            push_toml_number(
                &mut entries,
                "fragment_stride",
                args.fragment_stride,
                sources.fragment_stride,
            );
        }
        if sources.min_fragment_length.is_some() && args.min_fragment_length != args.fragment_length
        {
            push_toml_number(
                &mut entries,
                "min_fragment_length",
                args.min_fragment_length,
                sources.min_fragment_length,
            );
        }
        if sources.min_identity.is_some() && args.min_identity != DEFAULT_MIN_PERCENT_IDENTITY {
            push_toml_number(
                &mut entries,
                "min_identity",
                args.min_identity,
                sources.min_identity,
            );
        }
        if sources.mash_confidence.is_some() && args.mash_confidence != DEFAULT_MASH_CONFIDENCE {
            push_toml_number(
                &mut entries,
                "mash_confidence",
                args.mash_confidence,
                sources.mash_confidence,
            );
        }
        if sources.split_n_run.is_some() && args.split_n_run != DEFAULT_SPLIT_N_RUN {
            push_toml_number(
                &mut entries,
                "split_n_run",
                args.split_n_run,
                sources.split_n_run,
            );
        }
        if let Some(max_memory_bytes) = args.max_memory_bytes {
            push_toml_number(
                &mut entries,
                "max_memory_gb",
                max_memory_bytes as f64 / 1024.0 / 1024.0 / 1024.0,
                sources.max_memory_gb,
            );
        }
        if let Some(max_concurrent_shards) = args.max_concurrent_shards {
            push_toml_number(
                &mut entries,
                "max_concurrent_shards",
                max_concurrent_shards,
                sources.max_concurrent_shards,
            );
        }
        if sources.max_shard_minimizers.is_some()
            && args.max_shard_minimizers != DEFAULT_MAX_SHARD_MINIMIZERS
        {
            push_toml_number(
                &mut entries,
                "max_shard_minimizers",
                args.max_shard_minimizers,
                sources.max_shard_minimizers,
            );
        }
        if sources.index_build_mode.is_some() && args.index_build_mode != IndexBuildMode::Auto {
            push_toml_string(
                &mut entries,
                "index_build_mode",
                Some(&StartupValue::new(
                    args.index_build_mode.name(),
                    sources.index_build_mode.expect("checked above"),
                )),
            );
        }

        if entries.is_empty() {
            return;
        }

        eprintln!("################# FasterANI non-default runtime parameters #################");
        for entry in entries {
            eprintln!("{entry}");
        }
        eprintln!("############################################################################");
    }
}

#[derive(Default)]
struct ParameterSources {
    params_file: Option<ParameterSource>,
    threads: Option<ParameterSource>,
    freq_threshold_percent: Option<ParameterSource>,
    minmer_count: Option<ParameterSource>,
    kmer_size: Option<ParameterSource>,
    window_size: Option<ParameterSource>,
    fragment_length: Option<ParameterSource>,
    fragment_stride: Option<ParameterSource>,
    min_fragment_length: Option<ParameterSource>,
    min_identity: Option<ParameterSource>,
    mash_confidence: Option<ParameterSource>,
    split_n_run: Option<ParameterSource>,
    max_memory_gb: Option<ParameterSource>,
    max_shard_minimizers: Option<ParameterSource>,
    max_concurrent_shards: Option<ParameterSource>,
    index_build_mode: Option<ParameterSource>,
    reference_sketch: Option<ParameterSource>,
    tmp: Option<ParameterSource>,
    out: Option<ParameterSource>,
    mapping_stats: Option<ParameterSource>,
    bgzip: Option<ParameterSource>,
    header: Option<ParameterSource>,
    verbose: Option<ParameterSource>,
}

fn push_toml_array(entries: &mut Vec<String>, key: &str, values: &[StartupValue]) {
    if values.is_empty() {
        return;
    }

    let rendered_values: Vec<String> = values
        .iter()
        .map(|value| format!("\"{}\"", toml_escape(&value.value)))
        .collect();
    entries.push(format!(
        "{key} = [{}]{}",
        rendered_values.join(", "),
        source_comment_for_values(values)
    ));
}

fn push_toml_path(
    entries: &mut Vec<String>,
    key: &str,
    path: Option<&PathBuf>,
    source: Option<ParameterSource>,
) {
    let Some(path) = path else {
        return;
    };
    entries.push(format!(
        "{key} = \"{}\"{}",
        toml_escape(&path.display().to_string()),
        source_comment(source)
    ));
}

fn push_toml_string(entries: &mut Vec<String>, key: &str, value: Option<&StartupValue>) {
    let Some(value) = value else {
        return;
    };
    entries.push(format!(
        "{key} = \"{}\"{}",
        toml_escape(&value.value),
        source_comment(Some(value.source))
    ));
}

fn push_toml_bool(
    entries: &mut Vec<String>,
    key: &str,
    value: bool,
    source: Option<ParameterSource>,
) {
    if value {
        entries.push(format!("{key} = true{}", source_comment(source)));
    }
}

fn push_toml_number(
    entries: &mut Vec<String>,
    key: &str,
    value: impl std::fmt::Display,
    source: Option<ParameterSource>,
) {
    entries.push(format!("{key} = {value}{}", source_comment(source)));
}

fn source_comment(source: Option<ParameterSource>) -> String {
    source
        .map(|source| format!("  # from {}", source.label()))
        .unwrap_or_default()
}

fn source_comment_for_values(values: &[StartupValue]) -> String {
    let has_params_file = values
        .iter()
        .any(|value| value.source == ParameterSource::ParamsFile);
    let has_cli = values
        .iter()
        .any(|value| value.source == ParameterSource::Cli);

    match (has_params_file, has_cli) {
        (true, true) => "  # from params file + CLI".to_owned(),
        (true, false) => source_comment(Some(ParameterSource::ParamsFile)),
        (false, true) => source_comment(Some(ParameterSource::Cli)),
        (false, false) => String::new(),
    }
}

fn toml_escape(value: &str) -> String {
    let mut escaped = String::with_capacity(value.len());
    for ch in value.chars() {
        match ch {
            '"' => escaped.push_str("\\\""),
            '\\' => escaped.push_str("\\\\"),
            '\n' => escaped.push_str("\\n"),
            '\r' => escaped.push_str("\\r"),
            '\t' => escaped.push_str("\\t"),
            _ => escaped.push(ch),
        }
    }
    escaped
}

fn extract_params_file_path(args: &[String]) -> io::Result<Option<String>> {
    let mut params_file_path: Option<String> = None;
    let mut index = 0usize;
    while index < args.len() {
        if args[index] == "--params-file" {
            let value = args.get(index + 1).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidInput, "--params-file requires a path")
            })?;
            if value.starts_with("--") {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "--params-file requires a path",
                ));
            }
            if params_file_path.is_some() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "--params-file may only be supplied once",
                ));
            }
            params_file_path = Some(value.clone());
            index += 2;
        } else {
            index += 1;
        }
    }

    Ok(params_file_path)
}

fn params_file_base_dir(params_file_path: &Path) -> Option<PathBuf> {
    let absolute_path = path::absolute(params_file_path).ok()?;
    absolute_path.parent().map(Path::to_path_buf)
}

fn add_reference_file(
    value: &str,
    source: ParameterSource,
    base_dir: Option<&Path>,
    references: &mut Vec<FastaInput>,
    startup_output: &mut RuntimeStartupOutput,
) -> io::Result<()> {
    let reference_absolute_path = resolve_path_from_base(value, base_dir)?;
    match fs::metadata(&reference_absolute_path) {
        Ok(meta) if meta.is_file() => {}
        Ok(_) => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "Reference path is not a file: {}",
                    reference_absolute_path.display()
                ),
            ));
        }
        Err(err) => {
            return Err(io::Error::new(
                err.kind(),
                format!(
                    "Cannot access reference file {}: {}",
                    reference_absolute_path.display(),
                    err
                ),
            ));
        }
    }

    startup_output.reference_files.push(StartupValue::new(
        reference_absolute_path.to_string_lossy().into_owned(),
        source,
    ));
    let input_path = match source {
        ParameterSource::Cli => value.to_owned(),
        ParameterSource::ParamsFile => reference_absolute_path.display().to_string(),
    };
    references.push(FastaInput::from_path(input_path));
    Ok(())
}

fn add_reference_list(
    value: &str,
    source: ParameterSource,
    base_dir: Option<&Path>,
    references: &mut Vec<FastaInput>,
    startup_output: &mut RuntimeStartupOutput,
) -> io::Result<()> {
    let (absolute_path, validated_paths) = validate_and_read_path_list_from_base(value, base_dir)?;
    match fs::exists(&absolute_path) {
        Ok(true) => {}
        Ok(false) => eprintln!(
            "ERROR: Reference list {} does not exist.",
            absolute_path.display()
        ),
        Err(e) => eprintln!("ERROR: Error loading reference list: {}", e),
    }

    startup_output.reference_lists.push(StartupValue::new(
        absolute_path.to_string_lossy().into_owned(),
        source,
    ));
    references.extend(validated_paths.into_iter().map(FastaInput::from_path));
    Ok(())
}

fn add_query_file(
    value: &str,
    source: ParameterSource,
    base_dir: Option<&Path>,
    queries: &mut Vec<FastaInput>,
    startup_output: &mut RuntimeStartupOutput,
) -> io::Result<()> {
    if is_stdin_path(value) {
        startup_output
            .query_files
            .push(StartupValue::new(value, source));
        queries.push(FastaInput::from_stdin(None));
        return Ok(());
    }

    let query_absolute_path = resolve_path_from_base(value, base_dir)?;
    match fs::metadata(&query_absolute_path) {
        Ok(meta) if meta.is_file() => {}
        Ok(_) => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "Query path is not a file: {}",
                    query_absolute_path.display()
                ),
            ));
        }
        Err(err) => {
            return Err(io::Error::new(
                err.kind(),
                format!(
                    "Cannot access query file {}: {}",
                    query_absolute_path.display(),
                    err
                ),
            ));
        }
    }
    startup_output.query_files.push(StartupValue::new(
        query_absolute_path.to_string_lossy().into_owned(),
        source,
    ));
    let input_path = match source {
        ParameterSource::Cli => value.to_owned(),
        ParameterSource::ParamsFile => query_absolute_path.display().to_string(),
    };
    queries.push(FastaInput::from_path(input_path));
    Ok(())
}

fn add_query_list(
    value: &str,
    source: ParameterSource,
    base_dir: Option<&Path>,
    queries: &mut Vec<FastaInput>,
    startup_output: &mut RuntimeStartupOutput,
) -> io::Result<()> {
    let (absolute_path, validated_paths) = validate_and_read_path_list_from_base(value, base_dir)?;
    match fs::exists(&absolute_path) {
        Ok(true) => {}
        Ok(false) => eprintln!(
            "ERROR: Query list {} does not exist.",
            absolute_path.display()
        ),
        Err(e) => eprintln!("ERROR: Error loading query list: {}", e),
    }

    startup_output.query_lists.push(StartupValue::new(
        absolute_path.to_string_lossy().into_owned(),
        source,
    ));
    queries.extend(validated_paths.into_iter().map(FastaInput::from_path));
    Ok(())
}

fn max_memory_gb_to_bytes(value: f64, source: &str) -> io::Result<u64> {
    if !value.is_finite() || value <= 0.0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            describe_field_parsing_error(
                source,
                &value.to_string(),
                "must be a positive finite number",
            ),
        ));
    }

    Ok((value * 1024.0 * 1024.0 * 1024.0) as u64)
}

/// Parse command-line arguments.
pub(crate) fn parse_cli_args() -> io::Result<Option<CliArgs>> {
    let raw_args: Vec<String> = env::args().skip(1).collect();
    let params_file_path = extract_params_file_path(&raw_args)?;
    let params_file_config = match params_file_path.as_deref() {
        Some(path) => load_params_file(path)?,
        None => ParamsFileConfig::default(),
    };
    let params_file_base_dir = params_file_path
        .as_deref()
        .and_then(|path| params_file_base_dir(Path::new(path)));

    let mut references: Vec<FastaInput> = Vec::new();
    let mut queries: Vec<FastaInput> = Vec::new();
    let mut startup_output = RuntimeStartupOutput::default();
    let mut sources = ParameterSources::default();
    let mut stdin_query_name: Option<String> = None;
    let mut sketch_path: Option<PathBuf> = None;
    let mut tmp_dir: Option<PathBuf> = None;
    let mut out_path: Option<PathBuf> = None;
    let mut mapping_stats_path: Option<PathBuf> = None;
    let mut bgzip: bool = false;
    let mut emit_header: bool = false;
    let mut verbose: bool = false;
    let mut threads: usize = 1usize;
    let mut freq_threshold_percent: f64 = DEFAULT_FREQ_THRESHOLD_PERCENT;
    let mut minmer_count: Option<usize> = None;
    let mut kmer_size: usize = DEFAULT_KMER_SIZE;
    let mut window_size: usize = DEFAULT_WINDOW_SIZE;
    let mut fragment_length: u32 = DEFAULT_FRAGMENT_LENGTH;
    let mut fragment_stride: u32 = DEFAULT_FRAGMENT_STRIDE;
    let mut fragment_stride_was_set: bool = false;
    let mut min_fragment_length: u32 = DEFAULT_MIN_FRAGMENT_LENGTH;
    let mut min_fragment_length_was_set: bool = false;
    let mut min_identity: f64 = DEFAULT_MIN_PERCENT_IDENTITY;
    let mut mash_confidence: f64 = DEFAULT_MASH_CONFIDENCE;
    let mut split_n_run: usize = DEFAULT_SPLIT_N_RUN;
    let mut max_memory_bytes: Option<u64> = None;
    let mut max_concurrent_shards: Option<usize> = None;
    let mut max_shard_minimizers: Option<usize> = None;
    let mut index_build_mode: IndexBuildMode = IndexBuildMode::Auto;

    if let Some(path) = params_file_path.as_deref() {
        let absolute_path = path::absolute(path)?;
        startup_output.params_file = Some(StartupValue::new(
            absolute_path.to_string_lossy(),
            ParameterSource::Cli,
        ));
        sources.params_file = Some(ParameterSource::Cli);
    }

    let params_file_base_dir = params_file_base_dir.as_deref();
    if let Some(reference_files) = params_file_config.reference_files.as_ref() {
        for reference_file in reference_files {
            add_reference_file(
                reference_file,
                ParameterSource::ParamsFile,
                params_file_base_dir,
                &mut references,
                &mut startup_output,
            )?;
        }
    }
    if let Some(reference_lists) = params_file_config.reference_lists.as_ref() {
        for reference_list in reference_lists {
            add_reference_list(
                reference_list,
                ParameterSource::ParamsFile,
                params_file_base_dir,
                &mut references,
                &mut startup_output,
            )?;
        }
    }
    if let Some(query_files) = params_file_config.query_files.as_ref() {
        for query_file in query_files {
            add_query_file(
                query_file,
                ParameterSource::ParamsFile,
                params_file_base_dir,
                &mut queries,
                &mut startup_output,
            )?;
        }
    }
    if let Some(query_lists) = params_file_config.query_lists.as_ref() {
        for query_list in query_lists {
            add_query_list(
                query_list,
                ParameterSource::ParamsFile,
                params_file_base_dir,
                &mut queries,
                &mut startup_output,
            )?;
        }
    }
    if let Some(query_name) = params_file_config.query_name.as_ref() {
        startup_output.query_name =
            Some(StartupValue::new(query_name, ParameterSource::ParamsFile));
        stdin_query_name = Some(query_name.clone());
    }
    if let Some(reference_sketch) = params_file_config.reference_sketch.as_ref() {
        sketch_path = Some(resolve_path_from_base(
            reference_sketch,
            params_file_base_dir,
        )?);
        sources.reference_sketch = Some(ParameterSource::ParamsFile);
    }
    if let Some(tmp) = params_file_config.tmp.as_ref() {
        tmp_dir = Some(resolve_path_from_base(tmp, params_file_base_dir)?);
        sources.tmp = Some(ParameterSource::ParamsFile);
    }
    if let Some(out) = params_file_config.out.as_ref() {
        out_path = Some(resolve_path_from_base(out, params_file_base_dir)?);
        sources.out = Some(ParameterSource::ParamsFile);
    }
    if let Some(mapping_stats) = params_file_config.mapping_stats.as_ref() {
        mapping_stats_path = Some(resolve_path_from_base(mapping_stats, params_file_base_dir)?);
        sources.mapping_stats = Some(ParameterSource::ParamsFile);
    }
    if let Some(value) = params_file_config.bgzip {
        bgzip = value;
        sources.bgzip = Some(ParameterSource::ParamsFile);
    }
    if let Some(value) = params_file_config.header {
        emit_header = value;
        sources.header = Some(ParameterSource::ParamsFile);
    }
    if let Some(value) = params_file_config.verbose {
        verbose = value;
        sources.verbose = Some(ParameterSource::ParamsFile);
    }
    if let Some(value) = params_file_config.threads {
        threads = value;
        sources.threads = Some(ParameterSource::ParamsFile);
    }
    if let Some(value) = params_file_config.freq_threshold_percent {
        freq_threshold_percent = value;
        sources.freq_threshold_percent = Some(ParameterSource::ParamsFile);
    }
    if let Some(value) = params_file_config.minmer_count {
        minmer_count = Some(value);
        sources.minmer_count = Some(ParameterSource::ParamsFile);
    }
    if let Some(value) = params_file_config.kmer_size {
        kmer_size = value;
        sources.kmer_size = Some(ParameterSource::ParamsFile);
    }
    if let Some(value) = params_file_config.window_size {
        window_size = value;
        sources.window_size = Some(ParameterSource::ParamsFile);
    }
    if let Some(value) = params_file_config.fragment_length {
        fragment_length = value;
        sources.fragment_length = Some(ParameterSource::ParamsFile);
    }
    if let Some(value) = params_file_config.fragment_stride {
        fragment_stride = value;
        fragment_stride_was_set = true;
        sources.fragment_stride = Some(ParameterSource::ParamsFile);
    }
    if let Some(value) = params_file_config.min_fragment_length {
        min_fragment_length = value;
        min_fragment_length_was_set = true;
        sources.min_fragment_length = Some(ParameterSource::ParamsFile);
    }
    if let Some(value) = params_file_config.min_identity {
        min_identity = value;
        sources.min_identity = Some(ParameterSource::ParamsFile);
    }
    if let Some(value) = params_file_config.mash_confidence {
        mash_confidence = value;
        sources.mash_confidence = Some(ParameterSource::ParamsFile);
    }
    if let Some(value) = params_file_config.split_n_run {
        split_n_run = value;
        sources.split_n_run = Some(ParameterSource::ParamsFile);
    }
    if let Some(value) = params_file_config.max_memory_gb {
        max_memory_bytes = Some(max_memory_gb_to_bytes(value, "max_memory_gb")?);
        sources.max_memory_gb = Some(ParameterSource::ParamsFile);
    }
    if let Some(value) = params_file_config.max_shard_minimizers {
        max_shard_minimizers = Some(value);
        sources.max_shard_minimizers = Some(ParameterSource::ParamsFile);
    }
    if let Some(value) = params_file_config.max_concurrent_shards {
        max_concurrent_shards = Some(value);
        sources.max_concurrent_shards = Some(ParameterSource::ParamsFile);
    }
    if let Some(value) = params_file_config.index_build_mode.as_ref() {
        index_build_mode = value.parse::<IndexBuildMode>()?;
        sources.index_build_mode = Some(ParameterSource::ParamsFile);
    }

    let mut args = raw_args.into_iter();

    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--params-file" => {
                let _ = args.next().ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidInput, "--params-file requires a path")
                })?;
            }
            "--reference" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidInput, "--reference requires a path")
                })?;
                add_reference_file(
                    &value,
                    ParameterSource::Cli,
                    None,
                    &mut references,
                    &mut startup_output,
                )?;
            }
            "--reference-list" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--reference-list requires a path",
                    )
                })?;
                add_reference_list(
                    &value,
                    ParameterSource::Cli,
                    None,
                    &mut references,
                    &mut startup_output,
                )?;
            }
            "--query" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidInput, "--query requires a path")
                })?;
                if !is_stdin_path(&value) {
                    if stdin_query_name.is_some() {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidInput,
                            "--query-name may only be used with `--query -`",
                        ));
                    }
                }
                add_query_file(
                    &value,
                    ParameterSource::Cli,
                    None,
                    &mut queries,
                    &mut startup_output,
                )?;
            }
            "--query-name" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidInput, "--query-name requires a value")
                })?;
                if stdin_query_name.is_some() {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--query-name may only be supplied once",
                    ));
                }
                startup_output.query_name =
                    Some(StartupValue::new(value.clone(), ParameterSource::Cli));
                stdin_query_name = Some(value);
            }
            "--query-list" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidInput, "--query-list requires a path")
                })?;
                add_query_list(
                    &value,
                    ParameterSource::Cli,
                    None,
                    &mut queries,
                    &mut startup_output,
                )?;
            }
            "--reference-sketch" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--reference-sketch requires a path",
                    )
                })?;
                sketch_path = Some(PathBuf::from(value));
                sources.reference_sketch = Some(ParameterSource::Cli);
            }
            "--bgzip" => {
                bgzip = true;
                sources.bgzip = Some(ParameterSource::Cli);
            }
            "--header" => {
                emit_header = true;
                sources.header = Some(ParameterSource::Cli);
            }
            "--kmer-size" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidInput, "--kmer-size requires a value")
                })?;
                kmer_size = value.parse::<usize>().map_err(|err| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("invalid --kmer-size value {value:?}: {err}"),
                    )
                })?;
                sources.kmer_size = Some(ParameterSource::Cli);
                validate_kmer_size(kmer_size)?;
            }
            "--window-size" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--window-size requires a value",
                    )
                })?;
                window_size = value.parse::<usize>().map_err(|err| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("invalid --window-size value {value:?}: {err}"),
                    )
                })?;
                sources.window_size = Some(ParameterSource::Cli);
                validate_window_size(window_size)?;
            }
            "--fragment-length" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--fragment-length requires a value",
                    )
                })?;
                fragment_length = value.parse::<u32>().map_err(|err| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("invalid --fragment-length value {value:?}: {err}"),
                    )
                })?;
                sources.fragment_length = Some(ParameterSource::Cli);
                validate_fragment_length(fragment_length)?;
            }
            "--min-identity" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--min-identity requires a value",
                    )
                })?;
                min_identity = value.parse::<f64>().map_err(|err| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("invalid --min-identity value {value:?}: {err}"),
                    )
                })?;
                sources.min_identity = Some(ParameterSource::Cli);
                validate_min_identity(min_identity)?;
            }
            "--mash-confidence" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--mash-confidence requires a value",
                    )
                })?;
                mash_confidence = value.parse::<f64>().map_err(|err| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("invalid --mash-confidence value {value:?}: {err}"),
                    )
                })?;
                sources.mash_confidence = Some(ParameterSource::Cli);
                validate_mash_confidence(mash_confidence)?;
            }
            "--max-concurrent-shards" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--max-concurrent-shards requires a value",
                    )
                })?;
                let n = value.parse::<usize>().map_err(|err| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("invalid --max-concurrent-shards value {value:?}: {err}"),
                    )
                })?;
                if n == 0 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--max-concurrent-shards must be at least 1",
                    ));
                }
                max_concurrent_shards = Some(n);
                sources.max_concurrent_shards = Some(ParameterSource::Cli);
            }
            "--max-shard-minimizers" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--max-shard-minimizers requires a value",
                    )
                })?;
                let parsed_max_shard_minimizers: usize = value.parse::<usize>().map_err(|err| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("invalid --max-shard-minimizers value {value:?}: {err}"),
                    )
                })?;
                sources.max_shard_minimizers = Some(ParameterSource::Cli);
                validate_max_shard_minimizers(parsed_max_shard_minimizers)?;
                max_shard_minimizers = Some(parsed_max_shard_minimizers);
            }
            "--tmp" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidInput, "--tmp requires a directory")
                })?;
                tmp_dir = Some(PathBuf::from(value));
                sources.tmp = Some(ParameterSource::Cli);
            }
            "--index-build-mode" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--index-build-mode requires a value",
                    )
                })?;
                index_build_mode = value.parse::<IndexBuildMode>()?;
                sources.index_build_mode = Some(ParameterSource::Cli);
            }
            "--out" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidInput, "--out requires a path")
                })?;
                out_path = Some(PathBuf::from(value));
                sources.out = Some(ParameterSource::Cli);
            }
            "--mapping-stats" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--mapping-stats requires a path",
                    )
                })?;
                mapping_stats_path = Some(PathBuf::from(value));
                sources.mapping_stats = Some(ParameterSource::Cli);
            }
            "--threads" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidInput, "--threads requires a value")
                })?;
                threads = value.parse::<usize>().map_err(|err| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("invalid --threads value {value:?}: {err}"),
                    )
                })?;
                if threads == 0 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--threads must be at least 1",
                    ));
                }
                sources.threads = Some(ParameterSource::Cli);
            }
            "--max-memory-gb" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--max-memory-gb requires a value",
                    )
                })?;
                max_memory_bytes = Some(parse_max_memory_gb(&value)?);
                sources.max_memory_gb = Some(ParameterSource::Cli);
            }
            "--freq-threshold-percent" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--freq-threshold-percent requires a value",
                    )
                })?;
                freq_threshold_percent = value.parse::<f64>().map_err(|err| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("invalid --freq-threshold-percent value {value:?}: {err}"),
                    )
                })?;
                if !(0.0..=100.0).contains(&freq_threshold_percent) {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--freq-threshold-percent must be between 0 and 100",
                    ));
                }
                sources.freq_threshold_percent = Some(ParameterSource::Cli);
            }
            "--minmer-count" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--minmer-count requires a value",
                    )
                })?;
                let parsed_count: usize = value.parse::<usize>().map_err(|err| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("invalid --minmer-count value {value:?}: {err}"),
                    )
                })?;
                if parsed_count == 0 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--minmer-count must be at least 1",
                    ));
                }
                minmer_count = Some(parsed_count);
                sources.minmer_count = Some(ParameterSource::Cli);
            }
            "--fragment-stride" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--fragment-stride requires a value",
                    )
                })?;
                fragment_stride = value.parse::<u32>().map_err(|err| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("invalid --fragment-stride value {value:?}: {err}"),
                    )
                })?;
                fragment_stride_was_set = true;
                sources.fragment_stride = Some(ParameterSource::Cli);
                if fragment_stride == 0 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--fragment-stride must be at least 1",
                    ));
                }
            }
            "--min-fraglen" | "--min-fragment-length" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("{arg} requires a value"),
                    )
                })?;
                min_fragment_length = value.parse::<u32>().map_err(|err| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("invalid {arg} value {value:?}: {err}"),
                    )
                })?;
                min_fragment_length_was_set = true;
                sources.min_fragment_length = Some(ParameterSource::Cli);
                if min_fragment_length == 0 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("{arg} must be at least 1"),
                    ));
                }
            }
            "--split-N" | "--split-n" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("{arg} requires a value"),
                    )
                })?;
                split_n_run = value.parse::<usize>().map_err(|err| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("invalid {arg} value {value:?}: {err}"),
                    )
                })?;
                sources.split_n_run = Some(ParameterSource::Cli);
            }
            "--verbose" => {
                verbose = true;
                sources.verbose = Some(ParameterSource::Cli);
            }
            "--help" | "-h" | "--h" | "help" | "-?" => {
                eprintln!("{}", usage());
                return Ok(None);
            }
            "--version" | "v" => {
                eprintln!("fasterANI {}", env!("CARGO_PKG_VERSION"));
                process::exit(0);
            }
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("unknown argument {arg:?}\n{}", usage()),
                ));
            }
        }
    }

    if references.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("missing --reference\n{}", usage()),
        ));
    }

    if queries.is_empty() && sketch_path.is_none() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "missing --query (omit queries when using --sketch for build-only mode)\n{}",
                usage()
            ),
        ));
    }

    if stdin_query_name.is_some() {
        let label: String = stdin_query_name.take().expect("checked above");
        let Some(query) = queries.iter_mut().find(|query| is_stdin_path(&query.open)) else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "--query-name requires `--query -`",
            ));
        };
        query.label = label;
    }

    let stdin_reference_count: usize = references
        .iter()
        .filter(|reference| is_stdin_path(&reference.open))
        .count();
    let stdin_query_count: usize = queries
        .iter()
        .filter(|query| is_stdin_path(&query.open))
        .count();
    if stdin_reference_count > 1 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "only one `--reference -` may be supplied per run",
        ));
    }
    if stdin_query_count > 1 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "only one `--query -` may be supplied per run",
        ));
    }
    if stdin_reference_count > 0 && stdin_query_count > 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "reference and query cannot both be read from stdin in one run",
        ));
    }

    if threads == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--threads must be at least 1",
        ));
    }
    if !(0.0..=100.0).contains(&freq_threshold_percent) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--freq-threshold-percent must be between 0 and 100",
        ));
    }
    if minmer_count == Some(0) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--minmer-count must be at least 1",
        ));
    }
    if max_concurrent_shards == Some(0) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--max-concurrent-shards must be at least 1",
        ));
    }
    if !fragment_stride_was_set {
        fragment_stride = fragment_length;
    }
    if !min_fragment_length_was_set {
        min_fragment_length = fragment_length;
    }
    validate_kmer_size(kmer_size)?;
    validate_window_size(window_size)?;
    validate_fragment_length(fragment_length)?;
    validate_min_identity(min_identity)?;
    validate_mash_confidence(mash_confidence)?;
    if fragment_stride == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--fragment-stride must be at least 1",
        ));
    }
    if min_fragment_length == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--min-fraglen must be at least 1",
        ));
    }
    if fragment_stride > fragment_length {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("--fragment-stride must be <= --fragment-length ({fragment_length})"),
        ));
    }
    if min_fragment_length > fragment_length {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("--min-fraglen must be <= --fragment-length ({fragment_length})"),
        ));
    }

    let max_shard_minimizers: usize = max_shard_minimizers
        .unwrap_or_else(|| default_max_shard_minimizers_for_runtime(threads, max_memory_bytes));
    validate_max_shard_minimizers(max_shard_minimizers)?;

    let cli_args = CliArgs {
        references,
        queries,
        sketch_path,
        tmp_dir,
        out_path,
        mapping_stats_path,
        bgzip,
        emit_header,
        verbose,
        threads,
        freq_threshold_percent,
        minmer_count,
        kmer_size,
        window_size,
        fragment_length,
        fragment_stride,
        min_fragment_length,
        min_identity,
        mash_confidence,
        split_n_run,
        max_memory_bytes,
        max_concurrent_shards,
        max_shard_minimizers,
        index_build_mode,
    };
    startup_output.emit(&cli_args, &sources);

    Ok(Some(cli_args))
}

fn parse_max_memory_gb(value: &str) -> io::Result<u64> {
    let gb: f64 = value.parse::<f64>().map_err(|err| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("invalid --max-memory-gb value {value:?}: {err}"),
        )
    })?;
    if !gb.is_finite() || gb <= 0.0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--max-memory-gb must be a positive finite number",
        ));
    }

    Ok((gb * 1024.0 * 1024.0 * 1024.0) as u64)
}
