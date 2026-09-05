//! Command-line argument parsing and the `--help` text.

mod effective_config;
mod help;
mod input;

use std::{
    collections::HashSet,
    env, io, path,
    path::{Path, PathBuf},
    process,
};

use crate::ani::{
    constants::{
        DEFAULT_FRAGMENT_LENGTH, DEFAULT_FRAGMENT_STRIDE, DEFAULT_FREQ_THRESHOLD_PERCENT,
        DEFAULT_KMER_SIZE, DEFAULT_MASH_CONFIDENCE, DEFAULT_MAX_SHARD_SIZE_BYTES,
        DEFAULT_MINIMIZER_HASH_SEED, DEFAULT_MIN_FRAGMENT_LENGTH, DEFAULT_MIN_PERCENT_IDENTITY,
        DEFAULT_MPHF_GAMMA, DEFAULT_SPLIT_N_RUN, DEFAULT_WINDOW_SIZE,
    },
    io_util::{is_stdin_path, FastaInput},
    params_file::{load_params_file, ParamsFileConfig},
    sketch::serialize::{legacy_sketch_path, manifest_path},
    validation::{
        validate_fragment_length, validate_kmer_size, validate_mash_confidence,
        validate_mash_threshold, validate_max_shard_size_bytes, validate_mphf_gamma,
        validate_window_size,
    },
};
use effective_config::{ParameterSource, ParameterSources, RuntimeStartupOutput, StartupValue};
use help::usage;
use input::{
    add_query_file, add_query_list, add_reference_file, add_reference_list, params_file_base_dir,
    resolve_path_from_base,
};

/// Parsed command-line arguments.
pub(crate) struct CliArgs {
    pub(crate) references: Vec<FastaInput>,
    pub(crate) queries: Vec<FastaInput>,
    pub(crate) sketch_path: Option<PathBuf>,
    pub(crate) tmp_dir: Option<PathBuf>,
    pub(crate) out_path: Option<PathBuf>,
    pub(crate) mapping_stats_path: Option<PathBuf>,
    pub(crate) per_contig: bool,
    pub(crate) emit_header: bool,
    pub(crate) verbose: bool,
    pub(crate) quiet: bool,
    pub(crate) threads: usize,
    pub(crate) freq_threshold_percent: f64,
    pub(crate) minmer_count: Option<usize>,
    pub(crate) kmer_size: usize,
    pub(crate) window_size: usize,
    pub(crate) minimizer_hash_seed: u32,
    pub(crate) fragment_length: u32,
    pub(crate) fragment_stride: u32,
    pub(crate) min_fragment_length: u32,
    pub(crate) mash_threshold: f64,
    pub(crate) mash_confidence: f64,
    pub(crate) mphf_gamma: f64,
    /// Minimum run-length of ambiguous `N` bases that splits a contig; `0` disables.
    pub(crate) split_n_run: usize,
    pub(crate) max_shard_size_bytes: u64,
    pub(crate) shard_filter: Option<HashSet<usize>>,
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

fn parse_shard_filter(s: &str) -> io::Result<HashSet<usize>> {
    let mut indices = HashSet::new();
    for token in s.split(',') {
        let token = token.trim();
        if token.is_empty() {
            continue;
        }
        if let Some((start, end)) = token.split_once('-') {
            let start: usize = start.trim().parse::<usize>().map_err(|err| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("invalid --shards range start {start:?}: {err}"),
                )
            })?;
            let end: usize = end.trim().parse::<usize>().map_err(|err| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("invalid --shards range end {end:?}: {err}"),
                )
            })?;
            if start > end {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("invalid --shards range {start}-{end}: start must be <= end"),
                ));
            }
            indices.extend(start..=end);
        } else {
            let index: usize = token.parse::<usize>().map_err(|err| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("invalid --shards value {token:?}: {err}"),
                )
            })?;
            indices.insert(index);
        }
    }
    if indices.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--shards must specify at least one shard index",
        ));
    }
    Ok(indices)
}

fn parse_byte_size(value: &str) -> io::Result<u64> {
    let digit_count = value.bytes().take_while(u8::is_ascii_digit).count();
    let (number, suffix) = value.split_at(digit_count);
    if number.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("invalid --max-shard-size value {value:?}"),
        ));
    }

    let number = number.parse::<u64>().map_err(|err| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("invalid --max-shard-size value {value:?}: {err}"),
        )
    })?;
    let multiplier: u64 = match suffix {
        "" | "B" => 1,
        "KiB" => 1024,
        "MiB" => 1024 * 1024,
        "GiB" => 1024 * 1024 * 1024,
        _ => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "invalid --max-shard-size suffix in {value:?}; expected B, KiB, MiB, or GiB"
                ),
            ));
        }
    };
    let bytes = number.checked_mul(multiplier).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("--max-shard-size value {value:?} exceeds the u64 byte limit"),
        )
    })?;
    validate_max_shard_size_bytes(bytes)?;
    Ok(bytes)
}

/// Parse command-line arguments.
pub(crate) fn parse_cli_args() -> io::Result<Option<CliArgs>> {
    parse_cli_args_from(env::args().skip(1))
}

/// Parse command-line arguments supplied by a caller rather than the process environment.
pub(crate) fn parse_cli_args_from<I, S>(args: I) -> io::Result<Option<CliArgs>>
where
    I: IntoIterator<Item = S>,
    S: Into<String>,
{
    let raw_args: Vec<String> = args.into_iter().map(Into::into).collect();
    if raw_args.is_empty() {
        eprintln!("{}", usage());
        return Ok(None);
    }
    let params_file_path = extract_params_file_path(&raw_args)?;
    let skip_validation = raw_args.iter().any(|arg| arg == "--skip-validation");
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
    let mut cli_query_name_seen: bool = false;
    let mut sketch_path: Option<PathBuf> = None;
    let mut tmp_dir: Option<PathBuf> = None;
    let mut out_path: Option<PathBuf> = None;
    let mut mapping_stats_path: Option<PathBuf> = None;
    let mut per_contig: bool = false;
    let mut emit_header: bool = false;
    let mut verbose: bool = false;
    let mut quiet: bool = false;
    let mut threads: usize = 1usize;
    let mut freq_threshold_percent: f64 = DEFAULT_FREQ_THRESHOLD_PERCENT;
    let mut minmer_count: Option<usize> = None;
    let mut kmer_size: usize = DEFAULT_KMER_SIZE;
    let mut window_size: usize = DEFAULT_WINDOW_SIZE;
    let mut minimizer_hash_seed: u32 = DEFAULT_MINIMIZER_HASH_SEED;
    let mut fragment_length: u32 = DEFAULT_FRAGMENT_LENGTH;
    let mut fragment_stride: u32 = DEFAULT_FRAGMENT_STRIDE;
    let mut fragment_stride_was_set: bool = false;
    let mut min_fragment_length: u32 = DEFAULT_MIN_FRAGMENT_LENGTH;
    let mut min_fragment_length_was_set: bool = false;
    let mut mash_threshold: f64 = DEFAULT_MIN_PERCENT_IDENTITY;
    let mut mash_confidence: f64 = DEFAULT_MASH_CONFIDENCE;
    let mut mphf_gamma: f64 = DEFAULT_MPHF_GAMMA;
    let mut split_n_run: usize = DEFAULT_SPLIT_N_RUN;
    let mut max_shard_size_bytes: Option<u64> = None;
    let mut shard_filter: Option<HashSet<usize>> = None;

    if let Some(path) = params_file_path.as_deref() {
        let absolute_path = path::absolute(path)?;
        startup_output.params_file = Some(StartupValue::new(
            absolute_path.to_string_lossy(),
            ParameterSource::Cli,
        ));
    }

    let params_file_base_dir = params_file_base_dir.as_deref();
    if let Some(reference_files) = params_file_config.reference_files.as_ref() {
        for reference_file in reference_files {
            add_reference_file(
                reference_file,
                ParameterSource::ParamsFile,
                params_file_base_dir,
                skip_validation,
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
                skip_validation,
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
                skip_validation,
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
                skip_validation,
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
    if let Some(value) = params_file_config.header {
        emit_header = value;
        sources.header = Some(ParameterSource::ParamsFile);
    }
    if let Some(value) = params_file_config.per_contig {
        per_contig = value;
        sources.per_contig = Some(ParameterSource::ParamsFile);
    }
    if let Some(value) = params_file_config.verbose {
        verbose = value;
        sources.verbose = Some(ParameterSource::ParamsFile);
    }
    if let Some(value) = params_file_config.quiet {
        quiet = value;
        sources.quiet = Some(ParameterSource::ParamsFile);
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
    if let Some(value) = params_file_config.minimizer_hash_seed {
        minimizer_hash_seed = value;
        sources.minimizer_hash_seed = Some(ParameterSource::ParamsFile);
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
    if let Some(value) = params_file_config.mash_threshold {
        mash_threshold = value;
        sources.mash_threshold = Some(ParameterSource::ParamsFile);
    }
    if let Some(value) = params_file_config.mash_confidence {
        mash_confidence = value;
        sources.mash_confidence = Some(ParameterSource::ParamsFile);
    }
    if let Some(value) = params_file_config.mphf_gamma {
        mphf_gamma = value;
        sources.mphf_gamma = Some(ParameterSource::ParamsFile);
    }
    if let Some(value) = params_file_config.split_n_run {
        split_n_run = value;
        sources.split_n_run = Some(ParameterSource::ParamsFile);
    }
    if let Some(value) = params_file_config.max_shard_size.as_deref() {
        max_shard_size_bytes = Some(parse_byte_size(value)?);
        sources.max_shard_size_bytes = Some(ParameterSource::ParamsFile);
    }
    if let Some(value) = params_file_config.shards.as_ref() {
        shard_filter = Some(parse_shard_filter(value)?);
        sources.shards = Some(ParameterSource::ParamsFile);
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
                    skip_validation,
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
                    skip_validation,
                    &mut references,
                    &mut startup_output,
                )?;
            }
            "--query" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidInput, "--query requires a path")
                })?;
                add_query_file(
                    &value,
                    ParameterSource::Cli,
                    None,
                    skip_validation,
                    &mut queries,
                    &mut startup_output,
                )?;
            }
            "--query-name" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidInput, "--query-name requires a value")
                })?;
                if cli_query_name_seen {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--query-name may only be supplied once",
                    ));
                }
                cli_query_name_seen = true;
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
                    skip_validation,
                    &mut queries,
                    &mut startup_output,
                )?;
            }
            "--skip-validation" => {}
            "--per-contig" => {
                per_contig = true;
                sources.per_contig = Some(ParameterSource::Cli);
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
            "--minimizer-hash-seed" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--minimizer-hash-seed requires a value",
                    )
                })?;
                minimizer_hash_seed = value.parse::<u32>().map_err(|err| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("invalid --minimizer-hash-seed value {value:?}: {err}"),
                    )
                })?;
                sources.minimizer_hash_seed = Some(ParameterSource::Cli);
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
            "--mash-threshold" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--mash-threshold requires a value",
                    )
                })?;
                mash_threshold = value.parse::<f64>().map_err(|err| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("invalid --mash-threshold value {value:?}: {err}"),
                    )
                })?;
                sources.mash_threshold = Some(ParameterSource::Cli);
                validate_mash_threshold(mash_threshold)?;
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
            "--mphf-gamma" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidInput, "--mphf-gamma requires a value")
                })?;
                mphf_gamma = value.parse::<f64>().map_err(|err| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("invalid --mphf-gamma value {value:?}: {err}"),
                    )
                })?;
                sources.mphf_gamma = Some(ParameterSource::Cli);
                validate_mphf_gamma(mphf_gamma)?;
            }
            "--shards" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidInput, "--shards requires a value")
                })?;
                shard_filter = Some(parse_shard_filter(&value)?);
                sources.shards = Some(ParameterSource::Cli);
            }
            "--max-shard-size" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--max-shard-size requires a value",
                    )
                })?;
                let parsed_max_shard_size_bytes = parse_byte_size(&value)?;
                sources.max_shard_size_bytes = Some(ParameterSource::Cli);
                max_shard_size_bytes = Some(parsed_max_shard_size_bytes);
            }
            "--tmp" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidInput, "--tmp requires a directory")
                })?;
                tmp_dir = Some(PathBuf::from(value));
                sources.tmp = Some(ParameterSource::Cli);
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
            "--max-reference-frequency" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--max-reference-frequency requires a value",
                    )
                })?;
                freq_threshold_percent = value.parse::<f64>().map_err(|err| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("invalid --max-reference-frequency value {value:?}: {err}"),
                    )
                })?;
                if !(0.0..=100.0).contains(&freq_threshold_percent) {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--max-reference-frequency must be between 0 and 100",
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
            "--quiet" | "--silent" => {
                quiet = true;
                sources.quiet = Some(ParameterSource::Cli);
            }
            "--help" | "-h" | "--h" | "help" | "-?" => {
                eprintln!("{}", usage());
                return Ok(None);
            }
            "--version" | "-v" => {
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

    if references.is_empty() && sketch_path.is_none() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("missing --reference or --reference-sketch\n{}", usage()),
        ));
    }

    if queries.is_empty() && sketch_path.is_none() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "missing --query (omit queries when using --reference-sketch for build-only mode)\n{}",
                usage()
            ),
        ));
    }

    let existing_sketch_requested: bool = sketch_path.as_deref().is_some_and(|prefix| {
        manifest_path(prefix).exists() || legacy_sketch_path(prefix).is_some()
    });
    if queries.is_empty() && existing_sketch_requested {
        eprintln!(
            "WARNING\treference sketch already exists; no query was provided, leaving it unchanged"
        );
        return Ok(None);
    }

    let stdin_reference_count: usize = references
        .iter()
        .filter(|reference| is_stdin_path(&reference.input_path))
        .count();
    let stdin_query_count: usize = queries
        .iter()
        .filter(|query| is_stdin_path(&query.input_path))
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
    if let Some(label) = stdin_query_name.take() {
        let Some(query) = queries
            .iter_mut()
            .find(|query| is_stdin_path(&query.input_path))
        else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "--query-name requires `--query -`",
            ));
        };
        query.output_label = label;
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
            "--max-reference-frequency must be between 0 and 100",
        ));
    }
    if minmer_count == Some(0) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--minmer-count must be at least 1",
        ));
    }
    if shard_filter.is_some() && sketch_path.is_none() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--shards requires --reference-sketch",
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
    validate_mash_threshold(mash_threshold)?;
    validate_mash_confidence(mash_confidence)?;
    validate_mphf_gamma(mphf_gamma)?;
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

    let max_shard_size_bytes: u64 = max_shard_size_bytes.unwrap_or(DEFAULT_MAX_SHARD_SIZE_BYTES);
    validate_max_shard_size_bytes(max_shard_size_bytes)?;

    let cli_args = CliArgs {
        references,
        queries,
        sketch_path,
        tmp_dir,
        out_path,
        mapping_stats_path,
        per_contig,
        emit_header,
        verbose,
        quiet,
        threads,
        freq_threshold_percent,
        minmer_count,
        kmer_size,
        window_size,
        minimizer_hash_seed,
        fragment_length,
        fragment_stride,
        min_fragment_length,
        mash_threshold,
        mash_confidence,
        mphf_gamma,
        split_n_run,
        max_shard_size_bytes,
        shard_filter,
    };
    if !cli_args.quiet {
        startup_output.emit(&cli_args, &sources, skip_validation);
    }

    Ok(Some(cli_args))
}

#[cfg(test)]
mod tests {
    use super::{parse_byte_size, parse_cli_args_from, parse_shard_filter};

    #[test]
    fn byte_sizes_accept_raw_and_binary_units() {
        assert_eq!(parse_byte_size("512").unwrap(), 512);
        assert_eq!(parse_byte_size("512B").unwrap(), 512);
        assert_eq!(parse_byte_size("512KiB").unwrap(), 512 * 1024);
        assert_eq!(parse_byte_size("512MiB").unwrap(), 512 * 1024 * 1024);
        assert_eq!(parse_byte_size("4GiB").unwrap(), 4 * 1024 * 1024 * 1024);
    }

    #[test]
    fn byte_sizes_reject_zero_unknown_units_and_overflow() {
        assert!(parse_byte_size("0").is_err());
        assert!(parse_byte_size("1GB").is_err());
        assert!(parse_byte_size("1.5GiB").is_err());
        assert!(parse_byte_size("18446744073709551615GiB").is_err());
    }

    #[test]
    fn parse_cli_args_from_accepts_an_argument_iterator() {
        let args = parse_cli_args_from([
            "--reference",
            "reference.fna",
            "--query",
            "query.fna",
            "--skip-validation",
            "--threads",
            "3",
            "--quiet",
        ])
        .unwrap()
        .unwrap();

        assert_eq!(args.references[0].output_label, "reference.fna");
        assert_eq!(args.queries[0].output_label, "query.fna");
        assert_eq!(args.threads, 3);
        assert!(args.quiet);
    }

    #[test]
    fn parse_cli_args_from_rejects_invalid_values_without_a_process() {
        let error = parse_cli_args_from([
            "--reference",
            "reference.fna",
            "--query",
            "query.fna",
            "--skip-validation",
            "--threads",
            "0",
        ])
        .err()
        .expect("zero threads should be rejected");

        assert_eq!(error.kind(), std::io::ErrorKind::InvalidInput);
        assert!(error.to_string().contains("--threads must be at least 1"));
    }

    #[test]
    fn removed_bgzip_flag_is_rejected() {
        let error = parse_cli_args_from(["--bgzip"])
            .err()
            .expect("removed --bgzip flag should be rejected");

        assert_eq!(error.kind(), std::io::ErrorKind::InvalidInput);
        assert!(error.to_string().contains("unknown argument \"--bgzip\""));
    }

    #[test]
    fn removed_index_build_mode_flag_is_rejected() {
        let error = parse_cli_args_from(["--index-build-mode", "hash"])
            .err()
            .expect("removed --index-build-mode flag should be rejected");

        assert_eq!(error.kind(), std::io::ErrorKind::InvalidInput);
        assert!(error
            .to_string()
            .contains("unknown argument \"--index-build-mode\""));
    }

    #[test]
    fn removed_force_flag_is_rejected() {
        let error = parse_cli_args_from(["--force"])
            .err()
            .expect("removed --force flag should be rejected");

        assert_eq!(error.kind(), std::io::ErrorKind::InvalidInput);
        assert!(error.to_string().contains("unknown argument \"--force\""));
    }

    #[test]
    fn replaced_max_shard_minimizers_flag_is_rejected() {
        let error = parse_cli_args_from(["--max-shard-minimizers", "1000"])
            .err()
            .expect("replaced --max-shard-minimizers flag should be rejected");

        assert_eq!(error.kind(), std::io::ErrorKind::InvalidInput);
        assert!(error
            .to_string()
            .contains("unknown argument \"--max-shard-minimizers\""));
    }

    #[test]
    fn parse_shard_filter_handles_single_index() {
        let filter = parse_shard_filter("3").unwrap();
        assert_eq!(filter, [3].into());
    }

    #[test]
    fn parse_shard_filter_handles_comma_list() {
        let filter = parse_shard_filter("1,3,5").unwrap();
        assert_eq!(filter, [1, 3, 5].into());
    }

    #[test]
    fn parse_shard_filter_handles_range() {
        let filter = parse_shard_filter("5-8").unwrap();
        assert_eq!(filter, [5, 6, 7, 8].into());
    }

    #[test]
    fn parse_shard_filter_handles_mixed() {
        let filter = parse_shard_filter("1,10-11,12").unwrap();
        assert_eq!(filter, [1, 10, 11, 12].into());
    }

    #[test]
    fn parse_shard_filter_rejects_inverted_range() {
        assert!(parse_shard_filter("8-5").is_err());
    }

    #[test]
    fn parse_shard_filter_rejects_empty() {
        assert!(parse_shard_filter("").is_err());
    }

    #[test]
    fn parse_shard_filter_rejects_non_numeric() {
        assert!(parse_shard_filter("1,foo,3").is_err());
    }
}
