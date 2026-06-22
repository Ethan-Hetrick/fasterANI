//! Command-line argument parsing and the `--help` text.

use std::{env, fs, io, path::PathBuf};

use crate::ani::{
    default_shard_minimizers_for_runtime, is_stdin_path, validate_fragment_length,
    validate_kmer_size, validate_mash_confidence, validate_min_identity, validate_shard_minimizers,
    validate_shard_size, validate_window_size, FastaInput, IndexBuildMode, DEFAULT_FRAGMENT_LENGTH,
    DEFAULT_FRAGMENT_STRIDE, DEFAULT_FREQ_THRESHOLD_PERCENT, DEFAULT_KMER_SIZE,
    DEFAULT_MASH_CONFIDENCE, DEFAULT_MIN_FRAGMENT_LENGTH, DEFAULT_MIN_PERCENT_IDENTITY,
    DEFAULT_SHARD_SIZE, DEFAULT_SPLIT_N_RUN, DEFAULT_WINDOW_SIZE,
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
    pub(crate) shard_size: usize,
    pub(crate) shard_minimizers: usize,
    pub(crate) index_build_mode: IndexBuildMode,
}

fn usage() -> &'static str {
    "usage: fasterANI (--reference <ref.fa> | --reference-list <refs.txt>)... \
[(--query <query.fa> | --query-list <queries.txt>)...] [options]

Inputs:
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
  --minmer-count <n>            Keep only the n smallest-hash minimizers ('minmers')
                                  per query fragment (default: keep all).
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
  --bgzip                       Treat sketch sidecar inputs as bgzip-compressed.
  --shard-size <n>              References per shard (count; default 10000).
  --shard-minimizers <n>        Minimizer budget per shard (count; default: memory-aware).
  --index-build-mode <mode>     auto | hash | partitioned (default auto).

Resources:
  --threads <n>                 Worker threads, >= 1 (default 1).
  --max-memory-gb <gb>          Soft memory ceiling in GB (default: unlimited).
  --tmp <dir>                   Directory for temporary shard files.
  -h, --help                    Show this help."
}

fn validate_and_read_path_list(list_path: &str) -> io::Result<Vec<String>> {
    let contents = fs::read_to_string(list_path)?;
    let mut valid_paths = Vec::new();

    for line in contents.lines() {
        let path = line.trim();
        if path.is_empty() || path.starts_with('#') {
            continue;
        }

        match fs::metadata(path) {
            Ok(meta) if meta.is_file() => {
                valid_paths.push(path.to_owned());
            }
            Ok(_) => {
                return Err(io::Error::new(io::ErrorKind::InvalidInput, format!("Path in list is not a file: {}", path)));
            }
            Err(e) => {
                return Err(io::Error::new(e.kind(), format!("Cannot access path '{}' from list: {}\n{}", path, list_path, e)));
            }
        }
    }
    Ok(valid_paths)
}

/// Parse command-line arguments.
pub(crate) fn parse_cli_args() -> io::Result<Option<CliArgs>> {
    let mut references: Vec<FastaInput> = Vec::new();
    let mut queries: Vec<FastaInput> = Vec::new();
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
    let mut shard_size: usize = DEFAULT_SHARD_SIZE;
    let mut shard_minimizers: Option<usize> = None;
    let mut index_build_mode: IndexBuildMode = IndexBuildMode::Auto;
    let mut args: std::iter::Skip<std::env::Args> = env::args().skip(1);

    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--reference" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidInput, "--reference requires a path")
                })?;

                match fs::exists(&value) {
                    Ok(true) => println!(">>> Reference file: {}", value),
                    Ok(false) => println!("ERROR: Reference file {} does not exist.", value),
                    Err(e) => println!("ERROR: Error loading reference: {}", e),
                }

                references.push(FastaInput::from_path(value));
            }
            "--reference-list" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--reference-list requires a path",
                    )
                })?;
                match fs::exists(&value) {
                    Ok(true) => println!(">>> Reference list: {}", value),
                    Ok(false) => println!("ERROR: Reference list {} does not exist.", value),
                    Err(e) => println!("ERROR: Error loading reference list: {}", e),
                }
                let validated_paths = validate_and_read_path_list(&value)?;

                references.extend(
                    validated_paths
                        .into_iter()
                        .map(FastaInput::from_path),
                );
            }
            "--query" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidInput, "--query requires a path")
                })?;
                if is_stdin_path(&value) {
                    queries.push(FastaInput::from_stdin(None));
                } else {
                    if stdin_query_name.is_some() {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidInput,
                            "--query-name may only be used with `--query -`",
                        ));
                    }
                    match fs::exists(&value) {
                        Ok(true) => println!(">>> Query file: {}", value),
                        Ok(false) => println!("ERROR: Query file {} does not exist.", value),
                        Err(e) => println!("ERROR: Error loading query: {}", e),
                    }
                    queries.push(FastaInput::from_path(value));
                }
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
                stdin_query_name = Some(value);
            }
            "--query-list" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidInput, "--query-list requires a path")
                })?;
                match fs::exists(&value) {
                    Ok(true) => println!(">>> Query list: {}", value),
                    Ok(false) => println!("ERROR: Query list {} does not exist.", value),
                    Err(e) => println!("ERROR: Error loading query list: {}", e),
                }
                match fs::exists(&value) {
                    Ok(true) => println!(">>> Reference list: {}", value),
                    Ok(false) => println!("ERROR: Reference list {} does not exist.", value),
                    Err(e) => println!("ERROR: Error loading reference list: {}", e),
                }

                let validated_paths = validate_and_read_path_list(&value)?;

                queries.extend(
                    validated_paths
                        .into_iter()
                        .map(FastaInput::from_path),
                );
            }
            "--reference-sketch" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidInput, "--reference-sketch requires a path")
                })?;
                sketch_path = Some(PathBuf::from(value));
            }
            "--bgzip" => {
                bgzip = true;
            }
            "--header" => {
                emit_header = true;
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
                validate_mash_confidence(mash_confidence)?;
            }
            "--shard-size" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidInput, "--shard-size requires a value")
                })?;
                shard_size = value.parse::<usize>().map_err(|err| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("invalid --shard-size value {value:?}: {err}"),
                    )
                })?;
                validate_shard_size(shard_size)?;
            }
            "--shard-minimizers" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--shard-minimizers requires a value",
                    )
                })?;
                let parsed_shard_minimizers: usize = value.parse::<usize>().map_err(|err| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("invalid --shard-minimizers value {value:?}: {err}"),
                    )
                })?;
                validate_shard_minimizers(parsed_shard_minimizers)?;
                shard_minimizers = Some(parsed_shard_minimizers);
            }
            "--tmp" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidInput, "--tmp requires a directory")
                })?;
                tmp_dir = Some(PathBuf::from(value));
            }
            "--index-build-mode" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--index-build-mode requires a value",
                    )
                })?;
                index_build_mode = value.parse::<IndexBuildMode>()?;
            }
            "--out" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidInput, "--out requires a path")
                })?;
                out_path = Some(PathBuf::from(value));
            }
            "--mapping-stats" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--mapping-stats requires a path",
                    )
                })?;
                mapping_stats_path = Some(PathBuf::from(value));
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
            }
            "--max-memory-gb" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--max-memory-gb requires a value",
                    )
                })?;
                max_memory_bytes = Some(parse_max_memory_gb(&value)?);
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
            }
            "--verbose" => {
                verbose = true;
            }
            "--help" | "-h" => {
                eprintln!("{}", usage());
                return Ok(None);
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

    let shard_minimizers: usize = shard_minimizers
        .unwrap_or_else(|| default_shard_minimizers_for_runtime(threads, max_memory_bytes));
    validate_shard_minimizers(shard_minimizers)?;

    Ok(Some(CliArgs {
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
        shard_size,
        shard_minimizers,
        index_build_mode,
    }))
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
