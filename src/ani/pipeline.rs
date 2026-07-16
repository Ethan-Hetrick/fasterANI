//! Top-level `run` orchestration and output writing.

#[cfg(debug_assertions)]
use std::mem::size_of;
use std::{
    collections::HashSet,
    fs, io,
    io::{BufReader, BufWriter, Write},
    sync::{mpsc, Arc},
    thread,
    time::Instant,
};

use noodles::fasta;

use crate::ani::{
    compact_reciprocal_best_mappings, emit_progress, fastani_compatible_fragment_mode,
    final_ani_computation, is_no_usable_fragments_error, legacy_sketch_path, manifest_path,
    map_query_batch_to_reference_parallel, open_fasta_reader, parse_cli_args, peak_rss_kb,
    performance_metrics_enabled, shard_entry_path, AniComputation, AniSummary, CliArgs,
    ContigAniSummary, MappingExecutor, MappingMetrics, MappingOutput, MappingResult,
    MappingResultKey, QueryFile, QueryFragment, ReferenceContigName, ReferenceFile,
    ReferenceSketch, RuntimeOptions, ScratchFile, ShardedBuildOptions, SketchDatabase,
    SketchParams,
};
#[cfg(debug_assertions)]
use crate::ani::{
    memory_mib, ReferenceMemoryEstimate, SEED_HIT_HISTOGRAM_OVERFLOW_LABEL,
    SEED_HIT_HISTOGRAM_UPPER_BOUNDS,
};

#[derive(Clone, Copy)]
struct PairSummaryStats {
    pub(crate) emitted_pairs: usize,
    pub(crate) ani_min: f64,
    pub(crate) ani_max: f64,
    pub(crate) af_min: f64,
    pub(crate) af_max: f64,
}

impl Default for PairSummaryStats {
    fn default() -> Self {
        Self {
            emitted_pairs: 0,
            ani_min: f64::NAN,
            ani_max: f64::NAN,
            af_min: f64::NAN,
            af_max: f64::NAN,
        }
    }
}

impl PairSummaryStats {
    fn record_pair(&mut self, ani: f64, aligned_fraction: f64) {
        if self.emitted_pairs == 0 {
            self.ani_min = ani;
            self.ani_max = ani;
            self.af_min = aligned_fraction;
            self.af_max = aligned_fraction;
        } else {
            self.ani_min = self.ani_min.min(ani);
            self.ani_max = self.ani_max.max(ani);
            self.af_min = self.af_min.min(aligned_fraction);
            self.af_max = self.af_max.max(aligned_fraction);
        }
        self.emitted_pairs += 1;
    }

    fn merge(&mut self, other: Self) {
        if other.emitted_pairs == 0 {
            return;
        }
        if self.emitted_pairs == 0 {
            *self = other;
            return;
        }

        self.emitted_pairs += other.emitted_pairs;
        self.ani_min = self.ani_min.min(other.ani_min);
        self.ani_max = self.ani_max.max(other.ani_max);
        self.af_min = self.af_min.min(other.af_min);
        self.af_max = self.af_max.max(other.af_max);
    }
}

struct QueryOutputStats {
    pub(crate) summary_elapsed: std::time::Duration,
    pub(crate) pair_stats: PairSummaryStats,
}

struct RawQueryMappingStats {
    pub(crate) mapping_results: Vec<MappingResult>,
    pub(crate) metrics: MappingMetrics,
    #[cfg(debug_assertions)]
    pub(crate) max_mapping_result_bytes: usize,
}

struct RawQueryMappingBatch {
    pub(crate) queries: Vec<RawQueryMappingStats>,
    pub(crate) mapping_elapsed: std::time::Duration,
}

struct LoadedShard {
    pub(crate) shard_index: usize,
    pub(crate) shard_offset: usize,
    pub(crate) sketch: ReferenceSketch,
}

struct PreparedQuery {
    pub(crate) query_index: usize,
    pub(crate) query_path: String,
    pub(crate) query_start: Instant,
    pub(crate) scratch: ScratchFile,
    pub(crate) spool_bytes: u64,
    pub(crate) owned_bytes: usize,
    pub(crate) fragment_count: usize,
    pub(crate) minimizer_count: usize,
    pub(crate) seed_minimizer_count: usize,
}

struct CollectedQuery {
    pub(crate) query_index: usize,
    pub(crate) query_path: String,
    pub(crate) query_start: Instant,
    pub(crate) collection_elapsed: std::time::Duration,
    pub(crate) query_file: QueryFile,
    pub(crate) owned_bytes: usize,
    pub(crate) fragment_count: usize,
    pub(crate) minimizer_count: usize,
    pub(crate) seed_minimizer_count: usize,
}

enum QueryPreparation {
    Ready(CollectedQuery),
    Skipped {
        collection_elapsed: std::time::Duration,
    },
}

#[derive(Default)]
struct QueryAccounting {
    pub(crate) fragments: usize,
    pub(crate) minimizers: usize,
    pub(crate) seed_minimizers: usize,
    #[cfg(debug_assertions)]
    pub(crate) total_owned_bytes: usize,
    #[cfg(debug_assertions)]
    pub(crate) max_owned_bytes: usize,
}

impl QueryAccounting {
    fn record(&mut self, query: &CollectedQuery, collect_memory_metrics: bool) {
        self.fragments += query.fragment_count;
        self.minimizers += query.minimizer_count;
        self.seed_minimizers += query.seed_minimizer_count;

        #[cfg(debug_assertions)]
        if collect_memory_metrics {
            self.total_owned_bytes = self.total_owned_bytes.saturating_add(query.owned_bytes);
            self.max_owned_bytes = self.max_owned_bytes.max(query.owned_bytes);
        }

        #[cfg(not(debug_assertions))]
        let _ = collect_memory_metrics;
    }
}

const MAX_QUERY_BATCH_SPOOL_BYTES: u64 = 256 * 1024 * 1024;
const MAX_QUERY_BATCH_OWNED_BYTES: usize = 256 * 1024 * 1024;

fn can_admit_query_owned_bytes(
    query_count: usize,
    current_owned_bytes: usize,
    candidate_owned_bytes: usize,
) -> bool {
    query_count == 0
        || current_owned_bytes.saturating_add(candidate_owned_bytes) <= MAX_QUERY_BATCH_OWNED_BYTES
}

impl PreparedQuery {
    fn load(&self) -> io::Result<QueryFile> {
        let file = fs::File::open(&self.scratch.path).map_err(|error| {
            io::Error::new(
                error.kind(),
                format!(
                    "failed to open prepared query spool for '{}': {error}",
                    self.query_path
                ),
            )
        })?;
        QueryFile::read_from(&mut BufReader::new(file)).map_err(|error| {
            io::Error::new(
                error.kind(),
                format!(
                    "failed to load prepared query spool for '{}': {error}",
                    self.query_path
                ),
            )
        })
    }
}

fn prepare_query(
    args: &CliArgs,
    query_index: usize,
    completed_queries: &mut usize,
    total_start: Instant,
    progress_enabled: bool,
) -> io::Result<QueryPreparation> {
    let query = &args.queries[query_index];
    let query_path: String = query.label.clone();
    if progress_enabled {
        emit_progress(
            "query",
            &format!(
                "event=start\tquery_index={}\tqueries_started={}\tquery_total={}\tpath={query_path}",
                query_index + 1,
                query_index + 1,
                args.queries.len()
            ),
            total_start,
        );
    }

    let query_start: Instant = Instant::now();
    let mut reader_query: fasta::io::Reader<Box<dyn io::BufRead>> = open_fasta_reader(&query.open)?;
    let query_file: QueryFile = match QueryFile::collect(
        &mut reader_query,
        args.kmer_size,
        args.window_size,
        args.minimizer_hash_seed,
        args.minmer_count,
        args.fragment_length,
        args.fragment_stride,
        args.min_fragment_length,
        args.split_n_run,
        args.per_contig,
    ) {
        Ok(query_file) => query_file,
        Err(error) if args.queries.len() > 1 && is_no_usable_fragments_error(&error) => {
            let collection_elapsed: std::time::Duration = query_start.elapsed();
            if !args.quiet {
                eprintln!(
                    "WARNING\tevent=query_skipped\treason=no_usable_fragments\tpath={query_path}"
                );
            }
            *completed_queries += 1;
            if progress_enabled {
                emit_progress(
                    "query",
                    &format!(
                        "event=skipped\tquery_index={}\tquery_done={}\tquery_total={}\tpath={query_path}",
                        query_index + 1,
                        *completed_queries,
                        args.queries.len()
                    ),
                    total_start,
                );
            }
            return Ok(QueryPreparation::Skipped { collection_elapsed });
        }
        Err(error) => return Err(error),
    };
    let collection_elapsed: std::time::Duration = query_start.elapsed();
    let owned_bytes: usize = query_file.estimated_owned_bytes();

    Ok(QueryPreparation::Ready(CollectedQuery {
        query_index,
        query_path,
        query_start,
        collection_elapsed,
        fragment_count: query_file.fragments.len(),
        minimizer_count: query_file.total_minimizers(),
        seed_minimizer_count: query_file.total_seed_minimizers(),
        owned_bytes,
        query_file,
    }))
}

fn spool_query_file(
    query_file: &QueryFile,
    tmp_dir: Option<&std::path::Path>,
    query_index: usize,
) -> io::Result<(ScratchFile, u64)> {
    let (scratch, file) = ScratchFile::create(tmp_dir, &format!("prepared-query-{query_index}"))?;
    let mut writer = BufWriter::new(file);
    query_file.write_to(&mut writer)?;
    writer.flush()?;
    drop(writer);
    let spool_bytes: u64 = fs::metadata(&scratch.path)?.len();
    Ok((scratch, spool_bytes))
}

#[allow(clippy::too_many_arguments)]
fn collect_query_mapping_batch(
    mapping_executor: &MappingExecutor,
    reference_sketch: &ReferenceSketch,
    query_files: &[&QueryFile],
    kmer_size: usize,
    window_size: usize,
    min_identity: f64,
    mash_confidence: f64,
    frequency_threshold: usize,
    performance_metrics_enabled: bool,
) -> io::Result<RawQueryMappingBatch> {
    let mapping_start: Instant = Instant::now();
    let outputs: Vec<MappingOutput> = map_query_batch_to_reference_parallel(
        mapping_executor,
        reference_sketch,
        query_files,
        kmer_size,
        window_size,
        min_identity,
        mash_confidence,
        frequency_threshold,
        performance_metrics_enabled,
    )?;
    let mapping_elapsed: std::time::Duration = mapping_start.elapsed();
    let queries: Vec<RawQueryMappingStats> = outputs
        .into_iter()
        .map(|mapping_output| {
            let mapping_results: Vec<MappingResult> = mapping_output.results;
            #[cfg(debug_assertions)]
            let max_mapping_result_bytes: usize =
                mapping_results.capacity() * size_of::<MappingResult>();
            RawQueryMappingStats {
                mapping_results,
                metrics: mapping_output.metrics,
                #[cfg(debug_assertions)]
                max_mapping_result_bytes,
            }
        })
        .collect();

    Ok(RawQueryMappingBatch {
        queries,
        mapping_elapsed,
    })
}

fn write_results_header(output: &mut dyn Write, per_contig: bool) -> io::Result<()> {
    if per_contig {
        writeln!(
            output,
            "query_file\treference_file\tquery_contig\teligible_fragments\tshared_fragments\tshared_bases\tANI\tmedian_ANI\tstddev\tMAD\tci_95_upper\tci_95_lower\tF99\tF80"
        )
    } else {
        writeln!(
            output,
            "query_file\treference_file\tANI\tAF\ttotal_fragments\tmedian_ANI\tstddev\tMAD\tci_95_upper\tci_95_lower\tF99\tF80"
        )
    }
}

fn aggregate_values(
    summary: &AniSummary,
    query_mapped_length: u64,
    fragment_length: u32,
) -> (f64, f64, f64) {
    let total_fragment_equivalents: f64 = query_mapped_length as f64 / f64::from(fragment_length);
    let aligned_fraction: f64 = if total_fragment_equivalents > 0.0 {
        let shared_fragment_equivalents: f64 =
            summary.shared_bases as f64 / f64::from(fragment_length);
        shared_fragment_equivalents / total_fragment_equivalents
    } else {
        f64::NAN
    };
    let ani: f64 = if summary.shared_bases > 0 {
        summary.weighted_identity_sum / summary.shared_bases as f64
    } else {
        f64::NAN
    };

    (ani, aligned_fraction, total_fragment_equivalents)
}

fn contig_ani(summary: &AniSummary) -> f64 {
    if summary.shared_bases > 0 {
        summary.weighted_identity_sum / summary.shared_bases as f64
    } else {
        f64::NAN
    }
}

fn write_aggregate_summary_comments(
    reference_files: &[ReferenceFile],
    summaries: &[AniSummary],
    query_path: &str,
    query_mapped_length: u64,
    fragment_length: u32,
    output: &mut dyn Write,
) -> io::Result<()> {
    writeln!(
        output,
        "# aggregate_summary_header\tquery_file\treference_file\tANI\tAF\ttotal_fragments\tmedian_ANI\tstddev\tMAD\tci_95_upper\tci_95_lower\tF99\tF80"
    )?;

    for (reference_file, summary) in reference_files.iter().zip(summaries) {
        let (ani, aligned_fraction, total_fragment_equivalents) =
            aggregate_values(summary, query_mapped_length, fragment_length);
        let stats = summary.distribution_stats;
        writeln!(
            output,
            "# aggregate_summary\t{query_path}\t{}\t{ani:.3}\t{aligned_fraction:.3}\t{total_fragment_equivalents:.2}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}",
            reference_file.path,
            stats.median,
            stats.stddev,
            stats.mad,
            stats.ci_95_upper,
            stats.ci_95_lower,
            stats.f99,
            stats.f80,
        )?;
    }

    Ok(())
}

fn write_mapping_stats_header(output: &mut dyn Write) -> io::Result<()> {
    writeln!(
        output,
        "query_file\treference_file\tquery_contig\treference_contig\tquery_fragment_id\tquery_start\tquery_end\treference_start\treference_end\tidentity\tquery_minimizer_count\treference_minimizer_count\tshared_minimizers\tunion_minimizers\tjaccard\tfragment_length\tis_reciprocal_best"
    )
}

// Kept with the hidden P99/P80 calculations so significance markers can be
// restored deliberately if those columns are reintroduced.
#[allow(dead_code)]
fn significance_stars(p_value: f64) -> &'static str {
    if p_value < 0.0001 {
        "****"
    } else if p_value < 0.001 {
        "***"
    } else if p_value < 0.01 {
        "**"
    } else if p_value < 0.05 {
        "*"
    } else {
        ""
    }
}

fn compare_mapping_stats_rows(left: &MappingResult, right: &MappingResult) -> std::cmp::Ordering {
    left.reference_file_id
        .cmp(&right.reference_file_id)
        .then_with(|| left.query_fragment_id.cmp(&right.query_fragment_id))
        .then_with(|| left.reference_contig_id.cmp(&right.reference_contig_id))
        .then_with(|| left.reference_start.cmp(&right.reference_start))
        .then_with(|| left.query_fragment_length.cmp(&right.query_fragment_length))
        .then_with(|| left.identity.total_cmp(&right.identity))
        .then_with(|| left.query_minimizer_count.cmp(&right.query_minimizer_count))
        .then_with(|| {
            left.reference_minimizer_count
                .cmp(&right.reference_minimizer_count)
        })
        .then_with(|| left.shared_minimizers.cmp(&right.shared_minimizers))
        .then_with(|| left.union_minimizers.cmp(&right.union_minimizers))
        .then_with(|| left.jaccard.total_cmp(&right.jaccard))
}

fn write_mapping_stats(
    reference_files: &[ReferenceFile],
    reference_contig_names: Option<&[ReferenceContigName]>,
    mapping_results: &[MappingResult],
    reciprocal_best_keys: &HashSet<MappingResultKey>,
    query_file: &QueryFile,
    query_path: &str,
    output: &mut dyn Write,
) -> io::Result<()> {
    let mut ordered_mappings: Vec<&MappingResult> = mapping_results.iter().collect();
    ordered_mappings.sort_by(|left, right| compare_mapping_stats_rows(left, right));

    for mapping in ordered_mappings {
        let reference_file: &ReferenceFile = reference_files
            .get(mapping.reference_file_id)
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "mapping references missing reference file id {}",
                        mapping.reference_file_id
                    ),
                )
            })?;
        let query_fragment: &QueryFragment = query_file
            .fragments
            .get(mapping.query_fragment_id)
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "mapping references missing query fragment id {}",
                        mapping.query_fragment_id
                    ),
                )
            })?;
        let query_contig: &str = query_file
            .contig_names
            .get(query_fragment.contig_id)
            .map_or("unknown", String::as_str);

        let reference_contig_opt: Option<&str> = reference_contig_names
            .and_then(|contigs| contigs.get(mapping.reference_contig_id))
            .map(|c| c.name.as_str());

        let reference_offset: u64 = reference_contig_names
            .and_then(|contigs| contigs.get(mapping.reference_contig_id))
            .map_or(0, |contig| u64::from(contig.segment_start));
        let reference_start: u64 =
            reference_offset.saturating_add(u64::from(mapping.reference_start));
        let reference_end: u64 =
            reference_start.saturating_add(u64::from(mapping.query_fragment_length));
        let is_reciprocal_best: bool =
            reciprocal_best_keys.contains(&MappingResultKey::from_mapping(mapping));

        write!(
            output,
            "{query_path}\t{}\t{query_contig}\t",
            reference_file.path
        )?;

        if let Some(name) = reference_contig_opt {
            write!(output, "{name}")?;
        } else {
            write!(output, "{}", mapping.reference_contig_id)?;
        }

        writeln!(
            output,
            "\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{}\t{}\t{}\t{}\t{:.8}\t{}\t{}",
            mapping.query_fragment_id,
            query_fragment.start,
            query_fragment.end,
            reference_start,
            reference_end,
            mapping.identity,
            mapping.query_minimizer_count,
            mapping.reference_minimizer_count,
            mapping.shared_minimizers,
            mapping.union_minimizers,
            mapping.jaccard,
            mapping.query_fragment_length,
            is_reciprocal_best
        )?;
    }

    Ok(())
}

fn write_per_contig_results(
    reference_files: &[ReferenceFile],
    contig_summaries: &[Vec<ContigAniSummary>],
    query_file: &QueryFile,
    query_path: &str,
    output: &mut dyn Write,
) -> io::Result<()> {
    for (reference_file, per_contig) in reference_files.iter().zip(contig_summaries) {
        for (contig_id, contig_name) in query_file.contig_names.iter().enumerate() {
            let contig_summary = per_contig
                .get(contig_id)
                .cloned()
                .unwrap_or_else(ContigAniSummary::default);
            let summary = contig_summary.summary;
            let stats = summary.distribution_stats;
            writeln!(
                output,
                "{query_path}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}",
                reference_file.path,
                contig_name,
                contig_summary.eligible_fragments,
                summary.shared_fragments,
                summary.shared_bases,
                contig_ani(&summary),
                stats.median,
                stats.stddev,
                stats.mad,
                stats.ci_95_upper,
                stats.ci_95_lower,
                stats.f99,
                stats.f80,
            )?;
        }
    }

    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn write_query_outputs(
    reference_files: &[ReferenceFile],
    reference_contig_names: Option<&[ReferenceContigName]>,
    mapping_results: Vec<MappingResult>,
    query_file: &QueryFile,
    query_path: &str,
    output: &mut dyn Write,
    mapping_stats_output: Option<&mut (dyn Write + '_)>,
    fragment_length: u32,
    per_contig: bool,
    emit_header: bool,
) -> io::Result<QueryOutputStats> {
    let summary_start: Instant = Instant::now();
    let raw_mapping_results: Option<Vec<MappingResult>> = mapping_stats_output
        .is_some()
        .then(|| mapping_results.clone());
    let ani_computation: AniComputation = final_ani_computation(
        mapping_results,
        query_file,
        reference_files.len(),
        fragment_length,
    );

    if let (Some(stats_output), Some(raw_mapping_results)) =
        (mapping_stats_output, raw_mapping_results.as_deref())
    {
        write_mapping_stats(
            reference_files,
            reference_contig_names,
            raw_mapping_results,
            &ani_computation.reciprocal_best_keys,
            query_file,
            query_path,
            stats_output,
        )?;
    }
    let summary_elapsed: std::time::Duration = summary_start.elapsed();

    let query_mapped_length: u64 = query_file.mapped_length();
    let mut pair_stats: PairSummaryStats = PairSummaryStats::default();

    if per_contig {
        write_aggregate_summary_comments(
            reference_files,
            &ani_computation.summaries,
            query_path,
            query_mapped_length,
            fragment_length,
            output,
        )?;
        if emit_header {
            write_results_header(output, per_contig)?;
        }
        write_per_contig_results(
            reference_files,
            &ani_computation.contig_summaries,
            query_file,
            query_path,
            output,
        )?;
    }

    if emit_header && !per_contig {
        write_results_header(output, per_contig)?;
    }

    for (reference_file, summary) in reference_files.iter().zip(ani_computation.summaries.iter()) {
        if summary.shared_fragments == 0 {
            continue;
        }

        let (ani, aligned_fraction, total_fragment_equivalents) =
            aggregate_values(summary, query_mapped_length, fragment_length);
        let stats = summary.distribution_stats;
        // P99/P80 are still computed in AniDistributionStats for future experimentation,
        // but are intentionally not reported while their interpretation is unsettled.
        if !per_contig {
            writeln!(
                output,
                "{query_path}\t{}\t{ani:.3}\t{aligned_fraction:.3}\t{total_fragment_equivalents:.2}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}",
                reference_file.path,
                stats.median,
                stats.stddev,
                stats.mad,
                stats.ci_95_upper,
                stats.ci_95_lower,
                stats.f99,
                stats.f80,
            )?;
        }
        pair_stats.record_pair(ani, aligned_fraction);
    }

    Ok(QueryOutputStats {
        summary_elapsed,
        pair_stats,
    })
}

/// Run the fasterANI command-line application.
pub fn run() -> io::Result<()> {
    run_started_at(Instant::now())
}

/// Run the fasterANI command-line application using a caller-provided process start time.
pub fn run_started_at(total_start: Instant) -> io::Result<()> {
    let Some(args) = parse_cli_args()? else {
        return Ok(());
    };
    let progress_enabled: bool = args.verbose;
    let performance_metrics_enabled: bool = performance_metrics_enabled(&args);
    let runtime_options: RuntimeOptions = RuntimeOptions::default()
        .with_progress_enabled(progress_enabled)
        .with_worker_threads(args.threads)
        .with_mphf_gamma(args.mphf_gamma);
    if !fastani_compatible_fragment_mode(
        args.fragment_length,
        args.fragment_stride,
        args.min_fragment_length,
    ) {
        eprintln!(
            "WARNING\tevent=adaptive_fragment_mode\tfragment_length={}\tfragment_stride={}\tmin_fragment_length={}",
            args.fragment_length, args.fragment_stride, args.min_fragment_length
        );
    }
    if progress_enabled {
        emit_progress(
            "parameters",
            &format!(
                "event=algorithm\tkmer_size={}\twindow_size={}\tminimizer_hash_seed={}\tfragment_length={}\tfragment_stride={}\tmin_fragment_length={}\tmash_threshold={:.6}\tmash_confidence={:.6}\tmphf_gamma={:.6}\tminmer_count={}\tfreq_threshold_percent={:.6}\tsplit_n_run={}",
                args.kmer_size,
                args.window_size,
                args.minimizer_hash_seed,
                args.fragment_length,
                args.fragment_stride,
                args.min_fragment_length,
                args.mash_threshold,
                args.mash_confidence,
                args.mphf_gamma,
                args.minmer_count.map_or_else(|| "disabled".to_owned(), |count| count.to_string()),
                args.freq_threshold_percent,
                args.split_n_run,
            ),
            total_start,
        );
    }

    let sketch_was_requested: bool = args.sketch_path.is_some();
    let existing_database_loaded: bool = args.sketch_path.as_deref().is_some_and(|prefix| {
        manifest_path(prefix).exists() || legacy_sketch_path(prefix).is_some()
    });

    let reference_start: Instant = Instant::now();
    let mapping_stats_requested: bool = args.mapping_stats_path.is_some();
    let reference_database: SketchDatabase = SketchDatabase::collect_or_load(
        &args.references,
        SketchParams {
            kmer_size: args.kmer_size,
            window_size: args.window_size,
            minimizer_hash_seed: args.minimizer_hash_seed,
            fragment_length: args.fragment_length,
            min_fragment_length: args.min_fragment_length,
            split_n_run: args.split_n_run,
        },
        args.sketch_path.as_deref(),
        ShardedBuildOptions {
            tmp_dir: args.tmp_dir.as_deref(),
            bgzip: args.bgzip,
            max_shard_minimizers: args.max_shard_minimizers,
            index_build_mode: args.index_build_mode,
            threads: args.threads,
            force_rebuild: args.force && args.queries.is_empty(),
        },
        mapping_stats_requested,
        runtime_options,
    )?;
    let reference_mode: &str =
        reference_database.mode_name(sketch_was_requested, existing_database_loaded);
    let reference_elapsed: std::time::Duration = reference_start.elapsed();
    let mapping_executor = MappingExecutor::new(args.threads)?;
    #[cfg(debug_assertions)]
    let reference_memory: Option<ReferenceMemoryEstimate> = match &reference_database {
        SketchDatabase::Single(reference_sketch) if performance_metrics_enabled => {
            Some(reference_sketch.memory_estimate())
        }
        _ => None,
    };
    #[cfg(debug_assertions)]
    let frequency_threshold_report: String = match &reference_database {
        SketchDatabase::Single(reference_sketch) => reference_sketch
            .index
            .frequency_threshold(args.freq_threshold_percent)
            .to_string(),
        SketchDatabase::Sharded {
            global_frequencies, ..
        } => global_frequencies
            .frequency_threshold(args.freq_threshold_percent)
            .to_string(),
    };

    let mut query_elapsed: std::time::Duration = std::time::Duration::ZERO;
    let mut query_spool_elapsed: std::time::Duration = std::time::Duration::ZERO;
    let mut mapping_elapsed: std::time::Duration = std::time::Duration::ZERO;
    let mut summary_elapsed: std::time::Duration = std::time::Duration::ZERO;
    let mut query_accounting: QueryAccounting = QueryAccounting::default();
    #[cfg(debug_assertions)]
    let mut max_mapping_result_bytes: usize = 0usize;
    let mut mapping_count: usize = 0usize;
    let mut mapping_detail_metrics: MappingMetrics = MappingMetrics::default();
    let mut emitted_pairs: usize = 0usize;
    let mut pair_summary_stats: PairSummaryStats = PairSummaryStats::default();
    let mut skipped_queries: usize = 0usize;
    let mut completed_queries: usize = 0usize;
    let mut output: Box<dyn Write> = match &args.out_path {
        Some(path) => Box::new(BufWriter::new(fs::File::create(path)?)),
        None => Box::new(BufWriter::new(io::stdout())),
    };
    let mut results_header_written: bool = false;
    let mut mapping_stats_output: Option<Box<dyn Write>> = match &args.mapping_stats_path {
        Some(path) => {
            let mut writer: Box<dyn Write> = Box::new(BufWriter::new(fs::File::create(path)?));
            write_mapping_stats_header(&mut *writer)?;
            Some(writer)
        }
        None => None,
    };

    if let SketchDatabase::Sharded {
        prefix,
        manifest,
        global_frequencies,
    } = &reference_database
    {
        let frequency_threshold: usize =
            global_frequencies.frequency_threshold(args.freq_threshold_percent);
        let active_shards: Vec<_> = match &args.shard_filter {
            Some(filter) => {
                let manifest_indices: HashSet<usize> = manifest
                    .shards
                    .iter()
                    .map(|shard| shard.shard_index)
                    .collect();
                let mut unknown: Vec<usize> = filter
                    .iter()
                    .copied()
                    .filter(|index| !manifest_indices.contains(index))
                    .collect();
                if !unknown.is_empty() {
                    unknown.sort_unstable();
                    let mut available: Vec<usize> = manifest_indices.iter().copied().collect();
                    available.sort_unstable();
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!(
                            "--shards specified unknown shard indices: {}; available: {}",
                            unknown
                                .iter()
                                .map(std::string::ToString::to_string)
                                .collect::<Vec<_>>()
                                .join(", "),
                            available
                                .iter()
                                .map(std::string::ToString::to_string)
                                .collect::<Vec<_>>()
                                .join(", ")
                        ),
                    ));
                }
                manifest
                    .shards
                    .iter()
                    .filter(|&shard| filter.contains(&shard.shard_index))
                    .cloned()
                    .collect()
            }
            None => manifest.shards.clone(),
        };

        let mut reference_file_offsets: Vec<usize> = Vec::with_capacity(active_shards.len());
        let mut next_reference_file_offset: usize = 0usize;
        let mut reference_contig_offsets: Vec<usize> = Vec::with_capacity(active_shards.len());
        let mut next_reference_contig_offset: usize = 0usize;
        for shard in &active_shards {
            reference_file_offsets.push(next_reference_file_offset);
            next_reference_file_offset += shard.reference_count;
            reference_contig_offsets.push(next_reference_contig_offset);
            next_reference_contig_offset += shard.reference_contigs;
        }

        let mut prepared_queries: Vec<PreparedQuery> = Vec::with_capacity(args.queries.len());
        for query_index in 0..args.queries.len() {
            match prepare_query(
                &args,
                query_index,
                &mut completed_queries,
                total_start,
                progress_enabled,
            )? {
                QueryPreparation::Skipped { collection_elapsed } => {
                    query_elapsed += collection_elapsed;
                    skipped_queries += 1;
                }
                QueryPreparation::Ready(query) => {
                    query_elapsed += query.collection_elapsed;
                    query_accounting.record(&query, performance_metrics_enabled);
                    let spool_start: Instant = Instant::now();
                    let (scratch, spool_bytes) = spool_query_file(
                        &query.query_file,
                        args.tmp_dir.as_deref(),
                        query.query_index,
                    )?;
                    query_spool_elapsed += spool_start.elapsed();
                    prepared_queries.push(PreparedQuery {
                        query_index: query.query_index,
                        query_path: query.query_path,
                        query_start: query.query_start,
                        scratch,
                        spool_bytes,
                        owned_bytes: query.owned_bytes,
                        fragment_count: query.fragment_count,
                        minimizer_count: query.minimizer_count,
                        seed_minimizer_count: query.seed_minimizer_count,
                    });
                }
            }
        }

        if !prepared_queries.is_empty() {
            let mut per_query_mapping_results: Vec<Vec<MappingResult>> =
                (0..prepared_queries.len()).map(|_| Vec::new()).collect();
            let mut per_query_mapping_metrics: Vec<MappingMetrics> = (0..prepared_queries.len())
                .map(|_| MappingMetrics::default())
                .collect();
            let mut all_reference_files: Vec<ReferenceFile> =
                Vec::with_capacity(next_reference_file_offset);
            let mut all_reference_contig_names: Option<Vec<ReferenceContigName>> =
                mapping_stats_requested.then(|| Vec::with_capacity(next_reference_contig_offset));

            let (tx, rx) = mpsc::sync_channel::<io::Result<LoadedShard>>(1);
            let loader_active_shards = active_shards.clone();
            let loader_prefix = prefix.clone();
            let loader_tmp_dir = args.tmp_dir.clone();
            let loader_global_frequencies = Arc::clone(global_frequencies);
            let loader_params = SketchParams {
                kmer_size: args.kmer_size,
                window_size: args.window_size,
                minimizer_hash_seed: args.minimizer_hash_seed,
                fragment_length: args.fragment_length,
                min_fragment_length: args.min_fragment_length,
                split_n_run: args.split_n_run,
            };
            let active_shard_count = active_shards.len();
            let loader_handle = thread::spawn(move || {
                for (shard_offset, shard) in loader_active_shards.iter().enumerate() {
                    if progress_enabled {
                        emit_progress(
                            "shard_load",
                            &format!(
                                "event=start\tshard={}\tshards_total={}\tfilename={}\tloader=prefetch",
                                shard.shard_index, active_shard_count, shard.filename
                            ),
                            total_start,
                        );
                    }
                    let load_result = ReferenceSketch::load(
                        &shard_entry_path(&loader_prefix, shard),
                        loader_params,
                        mapping_stats_requested,
                        loader_tmp_dir.as_deref(),
                        runtime_options,
                    )
                    .map(|mut sketch| {
                        sketch.global_frequencies = Some(Arc::clone(&loader_global_frequencies));
                        sketch.prefetch_sequential();
                        LoadedShard {
                            shard_index: shard.shard_index,
                            shard_offset,
                            sketch,
                        }
                    });

                    let stop_after_send: bool = load_result.is_err();
                    if tx.send(load_result).is_err() || stop_after_send {
                        break;
                    }
                }
            });

            let mut shards_done: usize = 0usize;
            let shard_processing_result: io::Result<()> = (|| {
                while let Ok(load_result) = rx.recv() {
                    let mut loaded: LoadedShard = load_result?;
                    let reference_file_offset: usize = reference_file_offsets[loaded.shard_offset];
                    let reference_contig_offset: usize =
                        reference_contig_offsets[loaded.shard_offset];
                    let mut shard_mapping_count: usize = 0usize;

                    let mut query_batch_start: usize = 0;
                    while query_batch_start < prepared_queries.len() {
                        let mut query_batch_end: usize = query_batch_start;
                        let mut batch_spool_bytes: u64 = 0;
                        let mut batch_owned_bytes: usize = 0;
                        while query_batch_end < prepared_queries.len()
                            && query_batch_end - query_batch_start < args.threads.max(1)
                        {
                            let candidate_bytes: u64 =
                                prepared_queries[query_batch_end].spool_bytes;
                            let candidate_owned_bytes: usize =
                                prepared_queries[query_batch_end].owned_bytes;
                            if query_batch_end > query_batch_start
                                && (batch_spool_bytes.saturating_add(candidate_bytes)
                                    > MAX_QUERY_BATCH_SPOOL_BYTES
                                    || !can_admit_query_owned_bytes(
                                        query_batch_end - query_batch_start,
                                        batch_owned_bytes,
                                        candidate_owned_bytes,
                                    ))
                            {
                                break;
                            }
                            batch_spool_bytes = batch_spool_bytes.saturating_add(candidate_bytes);
                            batch_owned_bytes =
                                batch_owned_bytes.saturating_add(candidate_owned_bytes);
                            query_batch_end += 1;
                        }

                        let spool_load_start: Instant = Instant::now();
                        let query_files: Vec<QueryFile> = prepared_queries
                            [query_batch_start..query_batch_end]
                            .iter()
                            .map(PreparedQuery::load)
                            .collect::<io::Result<_>>()?;
                        let spool_load_elapsed: std::time::Duration = spool_load_start.elapsed();
                        query_spool_elapsed += spool_load_elapsed;
                        let query_file_refs: Vec<&QueryFile> = query_files.iter().collect();
                        let raw_batch: RawQueryMappingBatch = collect_query_mapping_batch(
                            &mapping_executor,
                            &loaded.sketch,
                            &query_file_refs,
                            args.kmer_size,
                            args.window_size,
                            args.mash_threshold,
                            args.mash_confidence,
                            frequency_threshold,
                            performance_metrics_enabled,
                        )?;
                        mapping_elapsed += raw_batch.mapping_elapsed;

                        for (batch_offset, mut raw_stats) in
                            raw_batch.queries.into_iter().enumerate()
                        {
                            let query_slot: usize = query_batch_start + batch_offset;
                            let query_shard_mapping_count: usize = raw_stats.mapping_results.len();
                            for mapping in &mut raw_stats.mapping_results {
                                mapping.reference_file_id += reference_file_offset;
                                mapping.reference_contig_id += reference_contig_offset;
                            }
                            if !mapping_stats_requested {
                                raw_stats.mapping_results = compact_reciprocal_best_mappings(
                                    raw_stats.mapping_results,
                                    args.fragment_length,
                                );
                            }

                            mapping_count += query_shard_mapping_count;
                            shard_mapping_count += query_shard_mapping_count;
                            per_query_mapping_metrics[query_slot].merge(raw_stats.metrics.clone());
                            mapping_detail_metrics.merge(raw_stats.metrics);
                            #[cfg(debug_assertions)]
                            if performance_metrics_enabled {
                                max_mapping_result_bytes = max_mapping_result_bytes
                                    .max(raw_stats.max_mapping_result_bytes);
                            }
                            per_query_mapping_results[query_slot]
                                .append(&mut raw_stats.mapping_results);
                        }

                        query_batch_start = query_batch_end;
                    }

                    all_reference_files.append(&mut loaded.sketch.files);
                    if let (Some(ref mut master_contigs), Some(mut shard_contigs)) =
                        (&mut all_reference_contig_names, loaded.sketch.contig_names)
                    {
                        master_contigs.append(&mut shard_contigs);
                    }

                    shards_done += 1;
                    if progress_enabled {
                        emit_progress(
                            "shard_load",
                            &format!(
                                "event=complete\tshard={}\tshards_done={shards_done}\tshards_total={}\tmappings={shard_mapping_count}",
                                loaded.shard_index,
                                active_shards.len(),
                            ),
                            total_start,
                        );
                    }
                }

                Ok(())
            })();
            drop(rx);
            let loader_join_result: io::Result<()> = loader_handle
                .join()
                .map_err(|_| io::Error::other("shard loader thread panicked"));
            if let Err(error) = shard_processing_result {
                let _ = loader_join_result;
                return Err(error);
            }
            loader_join_result?;

            for (query_slot, prepared_query) in prepared_queries.iter().enumerate() {
                let spool_load_start: Instant = Instant::now();
                let query_file = prepared_query.load()?;
                query_spool_elapsed += spool_load_start.elapsed();
                let results: Vec<MappingResult> =
                    std::mem::take(&mut per_query_mapping_results[query_slot]);
                let stats: QueryOutputStats = write_query_outputs(
                    &all_reference_files,
                    all_reference_contig_names.as_deref(),
                    results,
                    &query_file,
                    &prepared_query.query_path,
                    &mut *output,
                    mapping_stats_output.as_deref_mut(),
                    args.fragment_length,
                    args.per_contig,
                    args.emit_header && !results_header_written,
                )?;
                if args.emit_header {
                    results_header_written = true;
                }
                summary_elapsed += stats.summary_elapsed;
                emitted_pairs += stats.pair_stats.emitted_pairs;
                pair_summary_stats.merge(stats.pair_stats);
                let query_candidate_count: usize =
                    per_query_mapping_metrics[query_slot].candidate_regions_scored;

                let query_total_time = prepared_query.query_start.elapsed();
                completed_queries += 1;
                if progress_enabled {
                    emit_progress(
                        "query",
                        &format!(
                            "event=complete\tquery_index={}\tquery_done={completed_queries}\tquery_total={}\tfragments={}\tminimizers={}\tseed_minimizers={}\tcandidates={query_candidate_count}\telapsed_ms={}",
                            prepared_query.query_index + 1,
                            args.queries.len(),
                            prepared_query.fragment_count,
                            prepared_query.minimizer_count,
                            prepared_query.seed_minimizer_count,
                            query_total_time.as_millis()
                        ),
                        total_start,
                    );
                }
            }
        }
    }

    if let SketchDatabase::Single(reference_sketch) = &reference_database {
        let frequency_threshold: usize = reference_sketch
            .index
            .frequency_threshold(args.freq_threshold_percent);
        let mut next_query_index: usize = 0;
        let mut pending_query: Option<CollectedQuery> = None;
        while next_query_index < args.queries.len() || pending_query.is_some() {
            let mut query_batch: Vec<CollectedQuery> = Vec::with_capacity(args.threads.max(1));
            let mut batch_owned_bytes: usize = 0;
            while (next_query_index < args.queries.len() || pending_query.is_some())
                && query_batch.len() < args.threads.max(1)
            {
                let query: CollectedQuery = if let Some(query) = pending_query.take() {
                    query
                } else {
                    let query_index: usize = next_query_index;
                    next_query_index += 1;
                    match prepare_query(
                        &args,
                        query_index,
                        &mut completed_queries,
                        total_start,
                        progress_enabled,
                    )? {
                        QueryPreparation::Skipped { collection_elapsed } => {
                            query_elapsed += collection_elapsed;
                            skipped_queries += 1;
                            continue;
                        }
                        QueryPreparation::Ready(query) => {
                            query_elapsed += query.collection_elapsed;
                            query_accounting.record(&query, performance_metrics_enabled);
                            query
                        }
                    }
                };
                let prospective_owned_bytes: usize =
                    batch_owned_bytes.saturating_add(query.owned_bytes);
                if !can_admit_query_owned_bytes(
                    query_batch.len(),
                    batch_owned_bytes,
                    query.owned_bytes,
                ) {
                    pending_query = Some(query);
                    break;
                }
                batch_owned_bytes = prospective_owned_bytes;
                query_batch.push(query);
                if batch_owned_bytes >= MAX_QUERY_BATCH_OWNED_BYTES {
                    break;
                }
            }
            if query_batch.is_empty() {
                continue;
            }

            let query_file_refs: Vec<&QueryFile> =
                query_batch.iter().map(|query| &query.query_file).collect();
            let raw_batch: RawQueryMappingBatch = collect_query_mapping_batch(
                &mapping_executor,
                reference_sketch,
                &query_file_refs,
                args.kmer_size,
                args.window_size,
                args.mash_threshold,
                args.mash_confidence,
                frequency_threshold,
                performance_metrics_enabled,
            )?;
            mapping_elapsed += raw_batch.mapping_elapsed;

            for (query, raw_stats) in query_batch.into_iter().zip(raw_batch.queries) {
                let query_mapping_count: usize = raw_stats.mapping_results.len();
                let query_candidate_count: usize = raw_stats.metrics.candidate_regions_scored;
                let output_stats: QueryOutputStats = write_query_outputs(
                    &reference_sketch.files,
                    reference_sketch.contig_names.as_deref(),
                    raw_stats.mapping_results,
                    &query.query_file,
                    &query.query_path,
                    &mut *output,
                    mapping_stats_output.as_deref_mut(),
                    args.fragment_length,
                    args.per_contig,
                    args.emit_header && !results_header_written,
                )?;
                if args.emit_header {
                    results_header_written = true;
                }
                summary_elapsed += output_stats.summary_elapsed;
                mapping_count += query_mapping_count;
                emitted_pairs += output_stats.pair_stats.emitted_pairs;
                pair_summary_stats.merge(output_stats.pair_stats);
                mapping_detail_metrics.merge(raw_stats.metrics);
                #[cfg(debug_assertions)]
                {
                    max_mapping_result_bytes =
                        max_mapping_result_bytes.max(raw_stats.max_mapping_result_bytes);
                }

                completed_queries += 1;
                if progress_enabled {
                    emit_progress(
                        "query",
                        &format!(
                            "event=complete\tquery_index={}\tquery_done={completed_queries}\tquery_total={}\tfragments={}\tminimizers={}\tseed_minimizers={}\tcandidates={query_candidate_count}\telapsed_ms={}",
                            query.query_index + 1,
                            args.queries.len(),
                            query.fragment_count,
                            query.minimizer_count,
                            query.seed_minimizer_count,
                            query.query_start.elapsed().as_millis()
                        ),
                        total_start,
                    );
                }
            }
        }
    }

    output.flush()?;
    if let Some(mut stats_out) = mapping_stats_output {
        stats_out.flush()?;
    }

    let processed_queries: usize = args.queries.len().saturating_sub(skipped_queries);
    if !args.quiet {
        let peak_rss_kb_value: i64 = peak_rss_kb();
        let peak_rss_gib: f64 = if peak_rss_kb_value > 0 {
            peak_rss_kb_value as f64 / 1024.0 / 1024.0
        } else {
            f64::NAN
        };
        eprintln!(
            "SUMMARY\tstage=execution\tevent=complete\treference_mode={reference_mode}\tthreads={}\tqueries_requested={}\tqueries_processed={processed_queries}\tqueries_skipped={skipped_queries}\treferences={}\ttotal_runtime_s={:.3}\tpeak_rss_gib={peak_rss_gib:.3}\temitted_pairs={emitted_pairs}\tcandidate_regions={}\tani_min={:.3}\tani_max={:.3}\taf_min={:.3}\taf_max={:.3}\treference_build_s={:.3}\tmapping_s={:.3}\tquery_collect_s={:.3}\tquery_spool_io_s={:.3}\tcandidate_discovery_s={:.3}",
            args.threads,
            args.queries.len(),
            reference_database.reference_count(),
            total_start.elapsed().as_secs_f64(),
            mapping_detail_metrics.candidate_regions_scored,
            pair_summary_stats.ani_min,
            pair_summary_stats.ani_max,
            pair_summary_stats.af_min,
            pair_summary_stats.af_max,
            reference_elapsed.as_secs_f64(),
            mapping_elapsed.as_secs_f64(),
            query_elapsed.as_secs_f64(),
            query_spool_elapsed.as_secs_f64(),
            mapping_detail_metrics
                .candidate_discovery_elapsed
                .as_secs_f64(),
        );
    }

    if progress_enabled {
        emit_progress(
            "complete",
            &format!(
                "event=run\treference_mode={reference_mode}\tthreads={}\treferences={}\tqueries_requested={}\tqueries_processed={processed_queries}\tqueries_skipped={skipped_queries}\treference_contigs={}\tunique_minimizers={}\tquery_fragments={}\tquery_minimizers={}\tseed_minimizers={}\tcandidates={}\tmappings={mapping_count}\temitted_pairs={emitted_pairs}",
                args.threads,
                reference_database.reference_count(),
                args.queries.len(),
                reference_database.contig_count(),
                reference_database.unique_minimizer_count(),
                query_accounting.fragments,
                query_accounting.minimizers,
                query_accounting.seed_minimizers,
                mapping_detail_metrics.candidate_regions_scored
            ),
            total_start,
        );
    }

    #[cfg(not(debug_assertions))]
    let _ = (
        reference_elapsed,
        query_elapsed,
        query_spool_elapsed,
        mapping_elapsed,
        summary_elapsed,
        performance_metrics_enabled,
    );

    #[cfg(debug_assertions)]
    {
        if performance_metrics_enabled {
            eprintln!(
                "METRICS\treference_mode={reference_mode}\tthreads={}\tfreq_threshold_percent={:.6}\tfreq_threshold={}\tminmer_count={}\tfragment_stride={}\tmin_fragment_length={}\tsplit_n_run={}\treferences={}\tqueries={}\tskipped_queries={skipped_queries}\treference_contigs={}\tunique_minimizers={}\tquery_fragments={}\tmappings={mapping_count}\temitted_pairs={emitted_pairs}",
                args.threads,
                args.freq_threshold_percent,
                frequency_threshold_report,
                args.minmer_count.map_or_else(|| "disabled".to_owned(), |count| count.to_string()),
                args.fragment_stride,
                args.min_fragment_length,
                args.split_n_run,
                reference_database.reference_count(),
                args.queries.len(),
                reference_database.contig_count(),
                reference_database.unique_minimizer_count(),
                query_accounting.fragments,
            );
            eprintln!(
                "METRICS\treference_ms={:.3}\tquery_ms={:.3}\tquery_spool_io_ms={:.3}\tmapping_ms={:.3}\tsummary_ms={:.3}\ttotal_ms={:.3}\tpeak_rss_kb={}",
                reference_elapsed.as_secs_f64() * 1000.0,
                query_elapsed.as_secs_f64() * 1000.0,
                query_spool_elapsed.as_secs_f64() * 1000.0,
                mapping_elapsed.as_secs_f64() * 1000.0,
                summary_elapsed.as_secs_f64() * 1000.0,
                total_start.elapsed().as_secs_f64() * 1000.0,
                peak_rss_kb(),
            );
            eprintln!(
                "MAPPING_METRICS\tcandidate_discovery_calls={}\tseed_hits_collected={}\tcandidate_regions_found={}\tcandidate_regions_scored={}\treference_minimizers_scanned={}\tscoring_window_steps={}\tretained_mappings={}\tcandidate_discovery_ms={:.3}\tscoring_ms={:.3}",
                mapping_detail_metrics.candidate_discovery_calls,
                mapping_detail_metrics.seed_hits_collected,
                mapping_detail_metrics.candidate_regions_found,
                mapping_detail_metrics.candidate_regions_scored,
                mapping_detail_metrics.reference_minimizers_scanned,
                mapping_detail_metrics.scoring_window_steps,
                mapping_detail_metrics.retained_mappings,
                mapping_detail_metrics.candidate_discovery_elapsed.as_secs_f64() * 1000.0,
                mapping_detail_metrics.scoring_elapsed.as_secs_f64() * 1000.0,
            );
            eprintln!(
                "SEED_HIT_DISTRIBUTION\tseed_lookups={}\tzero_hit_lookups={}\tskipped_by_frequency={}\tmax_hit_list={}",
                mapping_detail_metrics.seed_lookup_count,
                mapping_detail_metrics.seed_lookup_zero_hits,
                mapping_detail_metrics.seed_lookup_skipped_by_frequency,
                mapping_detail_metrics.seed_hit_list_max,
            );
            for (bin_index, (&lookup_count, &hit_sum)) in mapping_detail_metrics
                .seed_hit_list_bins
                .iter()
                .zip(mapping_detail_metrics.seed_hit_list_bin_hits.iter())
                .enumerate()
            {
                let upper_bound: String =
                    SEED_HIT_HISTOGRAM_UPPER_BOUNDS.get(bin_index).map_or_else(
                        || SEED_HIT_HISTOGRAM_OVERFLOW_LABEL.to_owned(),
                        std::string::ToString::to_string,
                    );
                eprintln!(
                    "SEED_HIT_BIN\tupper={upper_bound}\tseed_lookups={lookup_count}\tseed_hits={hit_sum}"
                );
            }
            if let Some(reference_memory) = &reference_memory {
                eprintln!(
                    "MEMORY_ESTIMATE\tscope=reference\treference_minimizers={}\treference_minimizer_vec_bytes={}\treference_minimizer_vec_mib={:.3}\tunique_index_keys={}\tseed_hits={}\tseed_hit_vec_bytes={}\tseed_hit_vec_mib={:.3}\thash_index_rough_bytes={}\thash_index_rough_mib={:.3}\tmmap_file_bytes={}\tmmap_file_mib={:.3}",
                    reference_memory.reference_minimizers,
                    reference_memory.reference_minimizer_vec_bytes,
                    memory_mib(reference_memory.reference_minimizer_vec_bytes),
                    reference_memory.unique_index_keys,
                    reference_memory.seed_hits,
                    reference_memory.seed_hit_vec_bytes,
                    memory_mib(reference_memory.seed_hit_vec_bytes),
                    reference_memory.hash_index_rough_bytes,
                    memory_mib(reference_memory.hash_index_rough_bytes),
                    reference_memory.mmap_file_bytes,
                    memory_mib(reference_memory.mmap_file_bytes),
                );
            }
            eprintln!(
                "MEMORY_ESTIMATE\tscope=queries\ttotal_owned_bytes={}\ttotal_owned_mib={:.3}\tmax_query_owned_bytes={}\tmax_query_owned_mib={:.3}",
                query_accounting.total_owned_bytes,
                memory_mib(query_accounting.total_owned_bytes),
                query_accounting.max_owned_bytes,
                memory_mib(query_accounting.max_owned_bytes),
            );
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{
        can_admit_query_owned_bytes, compare_mapping_stats_rows, MappingResult,
        MAX_QUERY_BATCH_OWNED_BYTES,
    };

    fn mapping(reference_file_id: usize, query_fragment_id: usize, jaccard: f64) -> MappingResult {
        MappingResult {
            reference_file_id,
            reference_contig_id: 2,
            query_fragment_id,
            query_fragment_length: 3_000,
            reference_start: 12,
            identity: 98.5,
            query_minimizer_count: 20,
            reference_minimizer_count: 21,
            shared_minimizers: 19,
            union_minimizers: 22,
            jaccard,
        }
    }

    #[test]
    fn mapping_stats_rows_have_a_complete_deterministic_order() {
        let mut mappings = [
            mapping(1, 0, 0.9),
            mapping(0, 1, 0.8),
            mapping(0, 0, 0.7),
            mapping(0, 0, 0.6),
        ];

        mappings.sort_by(compare_mapping_stats_rows);

        let order: Vec<(usize, usize, u64)> = mappings
            .iter()
            .map(|row| {
                (
                    row.reference_file_id,
                    row.query_fragment_id,
                    row.jaccard.to_bits(),
                )
            })
            .collect();
        assert_eq!(
            order,
            vec![
                (0, 0, 0.6f64.to_bits()),
                (0, 0, 0.7f64.to_bits()),
                (0, 1, 0.8f64.to_bits()),
                (1, 0, 0.9f64.to_bits()),
            ]
        );
    }

    #[test]
    fn query_batch_memory_budget_allows_one_oversized_query_but_defers_the_next() {
        let two_hundred_mib: usize = 200 * 1024 * 1024;
        assert!(can_admit_query_owned_bytes(
            0,
            0,
            MAX_QUERY_BATCH_OWNED_BYTES + 1
        ));
        assert!(!can_admit_query_owned_bytes(
            1,
            two_hundred_mib,
            two_hundred_mib
        ));
        assert!(can_admit_query_owned_bytes(
            1,
            MAX_QUERY_BATCH_OWNED_BYTES - 1,
            1
        ));
    }
}
