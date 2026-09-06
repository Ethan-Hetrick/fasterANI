//! Sharded reference loading and multi-query execution.

use std::{
    collections::HashSet,
    io::{self, Write},
    sync::{mpsc, Arc},
    thread,
    time::Instant,
};

use crate::ani::{
    cli::CliArgs,
    mapping::{compact_reciprocal_best_mappings, MappingExecutor},
    metrics::MappingMetrics,
    model::{
        query::{MappingResult, QueryFile},
        reference::{ReferenceContigName, ReferenceFile, ReferenceSketch, SketchParams},
    },
    runtime::{emit_progress, RuntimeOptions},
    sketch::{database::SketchDatabase, serialize::shard_entry_path},
};

use super::{
    can_admit_query_owned_bytes, collect_query_mapping_batch, prepare_query, spool_query_file,
    write_query_outputs, LoadedShard, PreparedQuery, QueryOutputStats, QueryPreparation,
    RawQueryMappingBatch, RunStatistics, MAX_QUERY_BATCH_SPOOL_BYTES,
};

#[allow(clippy::too_many_arguments)]
pub(super) fn execute_sharded_queries(
    args: &CliArgs,
    reference_database: &SketchDatabase,
    mapping_executor: &MappingExecutor,
    mapping_stats_requested: bool,
    performance_metrics_enabled: bool,
    runtime_options: RuntimeOptions,
    total_start: Instant,
    progress_enabled: bool,
    output: &mut dyn Write,
    mapping_stats_output: &mut Option<Box<dyn Write>>,
    results_header_written: &mut bool,
    run_statistics: &mut RunStatistics,
    mapping_count: &mut usize,
    emitted_pairs: &mut usize,
    skipped_queries: &mut usize,
    completed_queries: &mut usize,
) -> io::Result<()> {
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
                completed_queries,
                total_start,
                progress_enabled,
            )? {
                QueryPreparation::Skipped { collection_elapsed } => {
                    run_statistics.query_elapsed += collection_elapsed;
                    *skipped_queries += 1;
                }
                QueryPreparation::Ready(query) => {
                    run_statistics.query_elapsed += query.collection_elapsed;
                    run_statistics
                        .query_accounting
                        .record(&query, performance_metrics_enabled);
                    let spool_start: Instant = Instant::now();
                    let (scratch, spool_bytes) = spool_query_file(
                        &query.query_file,
                        args.tmp_dir.as_deref(),
                        query.query_index,
                    )?;
                    run_statistics.query_spool_elapsed += spool_start.elapsed();
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
                        run_statistics.query_spool_elapsed += spool_load_elapsed;
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
                        run_statistics.mapping_elapsed += raw_batch.mapping_elapsed;

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

                            *mapping_count += query_shard_mapping_count;
                            shard_mapping_count += query_shard_mapping_count;
                            per_query_mapping_metrics[query_slot].merge(raw_stats.metrics.clone());
                            run_statistics
                                .mapping_detail_metrics
                                .merge(raw_stats.metrics);
                            #[cfg(debug_assertions)]
                            if performance_metrics_enabled {
                                run_statistics.max_mapping_result_bytes = run_statistics
                                    .max_mapping_result_bytes
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
                run_statistics.query_spool_elapsed += spool_load_start.elapsed();
                let results: Vec<MappingResult> =
                    std::mem::take(&mut per_query_mapping_results[query_slot]);
                let stats: QueryOutputStats = write_query_outputs(
                    &all_reference_files,
                    all_reference_contig_names.as_deref(),
                    results,
                    &query_file,
                    &prepared_query.query_path,
                    output,
                    mapping_stats_output.as_deref_mut(),
                    args.fragment_length,
                    args.per_contig,
                    args.emit_header && !*results_header_written,
                )?;
                if args.emit_header {
                    *results_header_written = true;
                }
                run_statistics.summary_elapsed += stats.summary_elapsed;
                *emitted_pairs += stats.pair_stats.emitted_pairs;
                run_statistics.pair_summary_stats.merge(stats.pair_stats);
                let query_candidate_count: usize =
                    per_query_mapping_metrics[query_slot].candidate_regions_scored;

                let query_total_time = prepared_query.query_start.elapsed();
                *completed_queries += 1;
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
    Ok(())
}
