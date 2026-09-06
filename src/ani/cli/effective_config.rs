//! Effective runtime-configuration provenance and startup rendering.

use std::path::PathBuf;

use super::CliArgs;

#[derive(Default)]
pub(super) struct RuntimeStartupOutput {
    pub(super) params_file: Option<StartupValue>,
    pub(super) reference_files: Vec<StartupValue>,
    pub(super) reference_lists: Vec<StartupValue>,
    pub(super) query_files: Vec<StartupValue>,
    pub(super) query_lists: Vec<StartupValue>,
    pub(super) query_name: Option<StartupValue>,
}

#[derive(Clone, Copy, Eq, PartialEq)]
pub(super) enum ParameterSource {
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

pub(super) struct StartupValue {
    value: String,
    source: ParameterSource,
}

impl StartupValue {
    pub(super) fn new(value: impl Into<String>, source: ParameterSource) -> Self {
        Self {
            value: value.into(),
            source,
        }
    }
}

impl RuntimeStartupOutput {
    pub(super) fn emit(&self, args: &CliArgs, sources: &ParameterSources, skip_validation: bool) {
        let mut entries: Vec<String> = Vec::new();

        push_cli_metadata_string(&mut entries, "params_file", self.params_file.as_ref());
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
        push_toml_bool(&mut entries, "header", args.emit_header, sources.header);
        push_toml_bool(
            &mut entries,
            "per_contig",
            args.per_contig,
            sources.per_contig,
        );
        push_toml_bool(&mut entries, "verbose", args.verbose, sources.verbose);
        push_toml_bool(&mut entries, "quiet", args.quiet, sources.quiet);
        entries.push(format!(
            "# skip_validation = {skip_validation}  (CLI-only metadata)"
        ));

        push_toml_number(&mut entries, "threads", args.threads, sources.threads);
        push_toml_number(
            &mut entries,
            "freq_threshold_percent",
            args.freq_threshold_percent,
            sources.freq_threshold_percent,
        );
        if let Some(minmer_count) = args.minmer_count {
            push_toml_number(
                &mut entries,
                "minmer_count",
                minmer_count,
                sources.minmer_count,
            );
        }
        push_toml_number(&mut entries, "kmer_size", args.kmer_size, sources.kmer_size);
        push_toml_number(
            &mut entries,
            "window_size",
            args.window_size,
            sources.window_size,
        );
        push_toml_number(
            &mut entries,
            "minimizer_hash_seed",
            args.minimizer_hash_seed,
            sources.minimizer_hash_seed,
        );
        push_toml_number(
            &mut entries,
            "fragment_length",
            args.fragment_length,
            sources.fragment_length,
        );
        push_toml_number(
            &mut entries,
            "fragment_stride",
            args.fragment_stride,
            sources.fragment_stride,
        );
        push_toml_number(
            &mut entries,
            "min_fragment_length",
            args.min_fragment_length,
            sources.min_fragment_length,
        );
        push_toml_number(
            &mut entries,
            "mash_threshold",
            args.mash_threshold,
            sources.mash_threshold,
        );
        push_toml_number(
            &mut entries,
            "mash_confidence",
            args.mash_confidence,
            sources.mash_confidence,
        );
        push_toml_number(
            &mut entries,
            "mphf_gamma",
            args.mphf_gamma,
            sources.mphf_gamma,
        );
        push_toml_number(
            &mut entries,
            "split_n_run",
            args.split_n_run,
            sources.split_n_run,
        );
        entries.push(format!(
            "max_shard_size = \"{}B\"{}",
            args.max_shard_size_bytes,
            source_comment(sources.max_shard_size_bytes)
        ));
        if let Some(filter) = &args.shard_filter {
            let mut sorted: Vec<usize> = filter.iter().copied().collect();
            sorted.sort_unstable();
            let rendered: Vec<String> = sorted
                .iter()
                .map(std::string::ToString::to_string)
                .collect();
            entries.push(format!(
                "shards = \"{}\"{}",
                rendered.join(","),
                source_comment(sources.shards)
            ));
        }

        eprintln!("################## FasterANI effective runtime parameters ##################");
        for entry in entries {
            eprintln!("{entry}");
        }
        eprintln!("############################################################################");
    }
}

#[derive(Default)]
pub(super) struct ParameterSources {
    pub(super) threads: Option<ParameterSource>,
    pub(super) freq_threshold_percent: Option<ParameterSource>,
    pub(super) minmer_count: Option<ParameterSource>,
    pub(super) kmer_size: Option<ParameterSource>,
    pub(super) window_size: Option<ParameterSource>,
    pub(super) minimizer_hash_seed: Option<ParameterSource>,
    pub(super) fragment_length: Option<ParameterSource>,
    pub(super) fragment_stride: Option<ParameterSource>,
    pub(super) min_fragment_length: Option<ParameterSource>,
    pub(super) mash_threshold: Option<ParameterSource>,
    pub(super) mash_confidence: Option<ParameterSource>,
    pub(super) mphf_gamma: Option<ParameterSource>,
    pub(super) split_n_run: Option<ParameterSource>,
    pub(super) max_shard_size_bytes: Option<ParameterSource>,
    pub(super) shards: Option<ParameterSource>,
    pub(super) reference_sketch: Option<ParameterSource>,
    pub(super) tmp: Option<ParameterSource>,
    pub(super) out: Option<ParameterSource>,
    pub(super) mapping_stats: Option<ParameterSource>,
    pub(super) header: Option<ParameterSource>,
    pub(super) per_contig: Option<ParameterSource>,
    pub(super) verbose: Option<ParameterSource>,
    pub(super) quiet: Option<ParameterSource>,
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

fn push_cli_metadata_string(entries: &mut Vec<String>, key: &str, value: Option<&StartupValue>) {
    let Some(value) = value else {
        return;
    };
    entries.push(format!(
        "# {key} = \"{}\"  ({}; CLI-only metadata)",
        toml_escape(&value.value),
        value.source.label()
    ));
}

fn push_toml_bool(
    entries: &mut Vec<String>,
    key: &str,
    value: bool,
    source: Option<ParameterSource>,
) {
    entries.push(format!("{key} = {value}{}", source_comment(source)));
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
