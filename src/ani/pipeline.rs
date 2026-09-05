//! Top-level `run` orchestration.

mod query_batch;
mod reporting;
mod sharded;

use std::{
    fs, io,
    io::{BufWriter, Write},
    time::Instant,
};

use query_batch::{
    can_admit_query_owned_bytes, collect_query_mapping_batch, prepare_query, spool_query_file,
    CollectedQuery, PreparedQuery, QueryAccounting, QueryPreparation, RawQueryMappingBatch,
    MAX_QUERY_BATCH_OWNED_BYTES, MAX_QUERY_BATCH_SPOOL_BYTES,
};
#[cfg(test)]
use reporting::compare_mapping_stats_rows;
use reporting::{
    write_mapping_stats_header, write_query_outputs, PairSummaryStats, QueryOutputStats,
};
use sharded::execute_sharded_queries;

#[cfg(test)]
use crate::ani::model::query::MappingResult;
use crate::ani::{
    cli::parse_cli_args,
    mapping::MappingExecutor,
    metrics::MappingMetrics,
    minimizer::fastani_compatible_fragment_mode,
    model::{
        query::QueryFile,
        reference::{ReferenceSketch, ShardedBuildOptions, SketchParams},
    },
    runtime::{emit_progress, peak_rss_kb, performance_metrics_enabled, RuntimeOptions},
    sketch::{
        database::SketchDatabase,
        serialize::{legacy_sketch_path, manifest_path},
    },
};
#[cfg(debug_assertions)]
use crate::ani::{
    metrics::{SEED_HIT_HISTOGRAM_OVERFLOW_LABEL, SEED_HIT_HISTOGRAM_UPPER_BOUNDS},
    model::reference::ReferenceMemoryEstimate,
    runtime::memory_mib,
};

#[derive(Default)]
struct RunStatistics {
    query_elapsed: std::time::Duration,
    query_spool_elapsed: std::time::Duration,
    mapping_elapsed: std::time::Duration,
    summary_elapsed: std::time::Duration,
    query_accounting: QueryAccounting,
    #[cfg(debug_assertions)]
    max_mapping_result_bytes: usize,
    mapping_detail_metrics: MappingMetrics,
    pair_summary_stats: PairSummaryStats,
}

struct LoadedShard {
    pub(crate) shard_index: usize,
    pub(crate) shard_offset: usize,
    pub(crate) sketch: ReferenceSketch,
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

    let mut run_statistics = RunStatistics::default();
    let mut mapping_count: usize = 0usize;
    let mut emitted_pairs: usize = 0usize;
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

    execute_sharded_queries(
        &args,
        &reference_database,
        &mapping_executor,
        mapping_stats_requested,
        performance_metrics_enabled,
        runtime_options,
        total_start,
        progress_enabled,
        &mut *output,
        &mut mapping_stats_output,
        &mut results_header_written,
        &mut run_statistics,
        &mut mapping_count,
        &mut emitted_pairs,
        &mut skipped_queries,
        &mut completed_queries,
    )?;

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
                            run_statistics.query_elapsed += collection_elapsed;
                            skipped_queries += 1;
                            continue;
                        }
                        QueryPreparation::Ready(query) => {
                            run_statistics.query_elapsed += query.collection_elapsed;
                            run_statistics
                                .query_accounting
                                .record(&query, performance_metrics_enabled);
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
            run_statistics.mapping_elapsed += raw_batch.mapping_elapsed;

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
                run_statistics.summary_elapsed += output_stats.summary_elapsed;
                mapping_count += query_mapping_count;
                emitted_pairs += output_stats.pair_stats.emitted_pairs;
                run_statistics
                    .pair_summary_stats
                    .merge(output_stats.pair_stats);
                run_statistics
                    .mapping_detail_metrics
                    .merge(raw_stats.metrics);
                #[cfg(debug_assertions)]
                {
                    run_statistics.max_mapping_result_bytes = run_statistics
                        .max_mapping_result_bytes
                        .max(raw_stats.max_mapping_result_bytes);
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
            run_statistics.mapping_detail_metrics.candidate_regions_scored,
            run_statistics.pair_summary_stats.ani_min,
            run_statistics.pair_summary_stats.ani_max,
            run_statistics.pair_summary_stats.af_min,
            run_statistics.pair_summary_stats.af_max,
            reference_elapsed.as_secs_f64(),
            run_statistics.mapping_elapsed.as_secs_f64(),
            run_statistics.query_elapsed.as_secs_f64(),
            run_statistics.query_spool_elapsed.as_secs_f64(),
            run_statistics.mapping_detail_metrics
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
                run_statistics.query_accounting.fragments,
                run_statistics.query_accounting.minimizers,
                run_statistics.query_accounting.seed_minimizers,
                run_statistics.mapping_detail_metrics.candidate_regions_scored
            ),
            total_start,
        );
    }

    #[cfg(not(debug_assertions))]
    let _ = (
        reference_elapsed,
        run_statistics.query_elapsed,
        run_statistics.query_spool_elapsed,
        run_statistics.mapping_elapsed,
        run_statistics.summary_elapsed,
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
                run_statistics.query_accounting.fragments,
            );
            eprintln!(
                "METRICS\treference_ms={:.3}\tquery_ms={:.3}\tquery_spool_io_ms={:.3}\tmapping_ms={:.3}\tsummary_ms={:.3}\ttotal_ms={:.3}\tpeak_rss_kb={}",
                reference_elapsed.as_secs_f64() * 1000.0,
                run_statistics.query_elapsed.as_secs_f64() * 1000.0,
                run_statistics.query_spool_elapsed.as_secs_f64() * 1000.0,
                run_statistics.mapping_elapsed.as_secs_f64() * 1000.0,
                run_statistics.summary_elapsed.as_secs_f64() * 1000.0,
                total_start.elapsed().as_secs_f64() * 1000.0,
                peak_rss_kb(),
            );
            eprintln!(
                "MAPPING_METRICS\tcandidate_discovery_calls={}\tseed_hits_collected={}\tcandidate_regions_found={}\tcandidate_regions_scored={}\treference_minimizers_scanned={}\tscoring_window_steps={}\tretained_mappings={}\tcandidate_discovery_ms={:.3}\tscoring_ms={:.3}",
                run_statistics.mapping_detail_metrics.candidate_discovery_calls,
                run_statistics.mapping_detail_metrics.seed_hits_collected,
                run_statistics.mapping_detail_metrics.candidate_regions_found,
                run_statistics.mapping_detail_metrics.candidate_regions_scored,
                run_statistics.mapping_detail_metrics.reference_minimizers_scanned,
                run_statistics.mapping_detail_metrics.scoring_window_steps,
                run_statistics.mapping_detail_metrics.retained_mappings,
                run_statistics.mapping_detail_metrics.candidate_discovery_elapsed.as_secs_f64() * 1000.0,
                run_statistics.mapping_detail_metrics.scoring_elapsed.as_secs_f64() * 1000.0,
            );
            eprintln!(
                "SEED_HIT_DISTRIBUTION\tseed_lookups={}\tzero_hit_lookups={}\tskipped_by_frequency={}\tmax_hit_list={}",
                run_statistics.mapping_detail_metrics.seed_lookup_count,
                run_statistics.mapping_detail_metrics.seed_lookup_zero_hits,
                run_statistics.mapping_detail_metrics.seed_lookup_skipped_by_frequency,
                run_statistics.mapping_detail_metrics.seed_hit_list_max,
            );
            for (bin_index, (&lookup_count, &hit_sum)) in run_statistics
                .mapping_detail_metrics
                .seed_hit_list_bins
                .iter()
                .zip(
                    run_statistics
                        .mapping_detail_metrics
                        .seed_hit_list_bin_hits
                        .iter(),
                )
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
                run_statistics.query_accounting.total_owned_bytes,
                memory_mib(run_statistics.query_accounting.total_owned_bytes),
                run_statistics.query_accounting.max_owned_bytes,
                memory_mib(run_statistics.query_accounting.max_owned_bytes),
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
