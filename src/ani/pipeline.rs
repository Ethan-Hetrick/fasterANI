//! Top-level `run` orchestration and output writing.

#[cfg(debug_assertions)]
use std::mem::size_of;
use std::{
    collections::HashSet,
    fs, io,
    io::{BufWriter, Write},
    sync::mpsc,
    thread,
    time::Instant,
};

use noodles::fasta;

use crate::ani::{
    emit_progress, fastani_compatible_fragment_mode, final_ani_computation,
    is_no_usable_fragments_error, legacy_sketch_path, manifest_path,
    map_query_to_reference_parallel, open_fasta_reader, parse_cli_args,
    performance_metrics_enabled, shard_entry_path, AniComputation, CliArgs, MappingOutput,
    MappingResult, MappingResultKey, QueryFile, QueryFragment, ReferenceContigName, ReferenceFile,
    ReferenceSketch, RuntimeOptions, ShardedBuildOptions, SketchDatabase, SketchParams,
};
#[cfg(debug_assertions)]
use crate::ani::{
    memory_mib, peak_rss_kb, MappingMetrics, QueryMemoryEstimate, ReferenceMemoryEstimate,
    SEED_HIT_HISTOGRAM_OVERFLOW_LABEL, SEED_HIT_HISTOGRAM_UPPER_BOUNDS,
};

struct QueryMappingStats {
    pub(crate) mapping_elapsed: std::time::Duration,
    pub(crate) summary_elapsed: std::time::Duration,
    pub(crate) mapping_count: usize,
    pub(crate) emitted_pairs: usize,
    #[cfg(debug_assertions)]
    pub(crate) max_mapping_result_bytes: usize,
    #[cfg(debug_assertions)]
    pub(crate) metrics: MappingMetrics,
}

struct RawQueryMappingStats {
    pub(crate) mapping_results: Vec<MappingResult>,
    pub(crate) mapping_elapsed: std::time::Duration,
    #[cfg(debug_assertions)]
    pub(crate) max_mapping_result_bytes: usize,
    #[cfg(debug_assertions)]
    pub(crate) metrics: MappingMetrics,
}

struct LoadedShard {
    pub(crate) shard_index: usize,
    pub(crate) shard_offset: usize,
    pub(crate) sketch: ReferenceSketch,
}

struct PreloadedQuery {
    pub(crate) query_index: usize,
    pub(crate) query_path: String,
    pub(crate) query_start: Instant,
    pub(crate) query_file: QueryFile,
    pub(crate) fragment_count: usize,
    pub(crate) minimizer_count: usize,
    pub(crate) seed_minimizer_count: usize,
}

#[allow(clippy::too_many_arguments)]
fn collect_query_mappings(
    reference_sketch: &ReferenceSketch,
    query_file: &QueryFile,
    kmer_size: usize,
    window_size: usize,
    min_identity: f64,
    mash_confidence: f64,
    mapping_threads: usize,
    frequency_threshold: usize,
    performance_metrics_enabled: bool,
) -> io::Result<RawQueryMappingStats> {
    let mapping_start: Instant = Instant::now();
    let mapping_output: MappingOutput = map_query_to_reference_parallel(
        reference_sketch,
        query_file,
        kmer_size,
        window_size,
        min_identity,
        mash_confidence,
        mapping_threads,
        frequency_threshold,
        performance_metrics_enabled,
    )?;
    let mapping_elapsed: std::time::Duration = mapping_start.elapsed();
    let mapping_results: Vec<MappingResult> = mapping_output.results;
    #[cfg(debug_assertions)]
    let max_mapping_result_bytes: usize = mapping_results.capacity() * size_of::<MappingResult>();
    #[cfg(debug_assertions)]
    let metrics: MappingMetrics = mapping_output.metrics;

    Ok(RawQueryMappingStats {
        mapping_results,
        mapping_elapsed,
        #[cfg(debug_assertions)]
        max_mapping_result_bytes,
        #[cfg(debug_assertions)]
        metrics,
    })
}

fn write_results_header(output: &mut dyn Write) -> io::Result<()> {
    writeln!(
        output,
        "query_file\treference_file\tANI\tAF\ttotal_fragments\tmedian_ANI\tstddev\tci_95_upper\tci_95_lower\tP99\tP80"
    )
}

fn write_mapping_stats_header(output: &mut dyn Write) -> io::Result<()> {
    writeln!(
        output,
        "query_file\treference_file\tquery_contig\treference_contig\tquery_fragment_id\tquery_start\tquery_end\treference_start\treference_end\tidentity\tquery_minimizer_count\treference_minimizer_count\tshared_minimizers\tunion_minimizers\tjaccard\tfragment_length\tis_reciprocal_best"
    )
}

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

fn write_mapping_stats(
    reference_files: &[ReferenceFile],
    reference_contig_names: Option<&[ReferenceContigName]>,
    mapping_results: &[MappingResult],
    reciprocal_best_keys: &HashSet<MappingResultKey>,
    query_file: &QueryFile,
    query_path: &str,
    output: &mut dyn Write,
) -> io::Result<()> {
    for mapping in mapping_results {
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
            .map(String::as_str)
            .unwrap_or("unknown");

        let reference_contig_opt: Option<&str> = reference_contig_names
            .and_then(|contigs| contigs.get(mapping.reference_contig_id))
            .map(|c| c.name.as_str());

        let reference_offset: u64 = reference_contig_names
            .and_then(|contigs| contigs.get(mapping.reference_contig_id))
            .map(|contig| u64::from(contig.segment_start))
            .unwrap_or(0);
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
            write!(output, "{}", name)?;
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
) -> io::Result<(std::time::Duration, usize)> {
    if let Some(stats_output) = mapping_stats_output {
        let reciprocal_keys = final_ani_computation(
            mapping_results.clone(),
            reference_files.len(),
            fragment_length,
        )
        .reciprocal_best_keys;
        write_mapping_stats(
            reference_files,
            reference_contig_names,
            &mapping_results,
            &reciprocal_keys,
            query_file,
            query_path,
            stats_output,
        )?;
    }

    let summary_start: Instant = Instant::now();
    let ani_computation: AniComputation =
        final_ani_computation(mapping_results, reference_files.len(), fragment_length);
    let summary_elapsed: std::time::Duration = summary_start.elapsed();

    let query_mapped_length: u64 = query_file.mapped_length();
    let mut emitted_pairs: usize = 0usize;

    for (reference_file, summary) in reference_files.iter().zip(ani_computation.summaries) {
        if summary.shared_fragments == 0 {
            continue;
        }

        let ani: f64 = summary.weighted_identity_sum / summary.shared_bases as f64;
        let shared_fragment_equivalents: f64 = summary.shared_bases as f64 / fragment_length as f64;
        let total_fragment_equivalents: f64 = query_mapped_length as f64 / fragment_length as f64;
        let aligned_fraction: f64 = shared_fragment_equivalents / total_fragment_equivalents;
        let stats = summary.distribution_stats;
        writeln!(
            output,
            "{query_path}\t{}\t{ani:.3}\t{aligned_fraction:.3}\t{total_fragment_equivalents:.2}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3e}{}\t{:.3e}{}",
            reference_file.path,
            stats.median,
            stats.stddev,
            stats.ci_95_upper,
            stats.ci_95_lower,
            stats.p99,
            significance_stars(stats.p99),
            stats.p80,
            significance_stars(stats.p80),
        )?;
        emitted_pairs += 1;
    }

    Ok((summary_elapsed, emitted_pairs))
}

#[allow(clippy::too_many_arguments)]
fn map_query_against_reference_sketch(
    reference_sketch: &ReferenceSketch,
    query_file: &QueryFile,
    query_path: &str,
    args: &CliArgs,
    frequency_threshold: usize,
    performance_metrics_enabled: bool,
    output: &mut dyn Write,
    mapping_stats_output: Option<&mut (dyn Write + '_)>,
) -> io::Result<QueryMappingStats> {
    let raw_stats: RawQueryMappingStats = collect_query_mappings(
        reference_sketch,
        query_file,
        args.kmer_size,
        args.window_size,
        args.mash_threshold,
        args.mash_confidence,
        args.threads,
        frequency_threshold,
        performance_metrics_enabled,
    )?;
    let mapping_count: usize = raw_stats.mapping_results.len();
    let (summary_elapsed, emitted_pairs) = write_query_outputs(
        &reference_sketch.files,
        reference_sketch.contig_names.as_deref(),
        raw_stats.mapping_results,
        query_file,
        query_path,
        output,
        mapping_stats_output,
        args.fragment_length,
    )?;

    Ok(QueryMappingStats {
        mapping_elapsed: raw_stats.mapping_elapsed,
        summary_elapsed,
        mapping_count,
        emitted_pairs,
        #[cfg(debug_assertions)]
        max_mapping_result_bytes: raw_stats.max_mapping_result_bytes,
        #[cfg(debug_assertions)]
        metrics: raw_stats.metrics,
    })
}

/// Run the fasterANI command-line application.
pub fn run() -> io::Result<()> {
    let total_start: Instant = Instant::now();
    let Some(args) = parse_cli_args()? else {
        return Ok(());
    };
    let progress_enabled: bool = args.verbose;
    let performance_metrics_enabled: bool = performance_metrics_enabled(&args);
    let runtime_options: RuntimeOptions = RuntimeOptions {
        progress_enabled,
        worker_threads: args.threads,
    };
    if !fastani_compatible_fragment_mode(
        args.fragment_length,
        args.fragment_stride,
        args.min_fragment_length,
    ) {
        eprintln!(
            "WARNING\tadaptive fragment mode enabled\tfragment_length={}\tfragment_stride={}\tmin_fragment_length={}",
            args.fragment_length, args.fragment_stride, args.min_fragment_length
        );
    }
    if progress_enabled {
        emit_progress(
            "parameters",
            &format!(
                "event=algorithm\tkmer_size={}\twindow_size={}\tfragment_length={}\tfragment_stride={}\tmin_fragment_length={}\tmash_threshold={:.6}\tmash_confidence={:.6}\tminmer_count={}\tfreq_threshold_percent={:.6}\tsplit_n_run={}",
                args.kmer_size,
                args.window_size,
                args.fragment_length,
                args.fragment_stride,
                args.min_fragment_length,
                args.mash_threshold,
                args.mash_confidence,
                args.minmer_count
                    .map(|count| count.to_string())
                    .unwrap_or_else(|| "disabled".to_owned()),
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
        },
        mapping_stats_requested,
        runtime_options,
    )?;
    let reference_mode: &str =
        reference_database.mode_name(sketch_was_requested, existing_database_loaded);
    let reference_elapsed: std::time::Duration = reference_start.elapsed();
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
        SketchDatabase::Sharded { .. } => "per-shard".to_string(),
    };

    let mut query_elapsed: std::time::Duration = std::time::Duration::ZERO;
    let mut mapping_elapsed: std::time::Duration = std::time::Duration::ZERO;
    let mut summary_elapsed: std::time::Duration = std::time::Duration::ZERO;
    let mut total_query_fragments_all: usize = 0usize;
    let mut total_query_minimizers_all: usize = 0usize;
    let mut total_query_seed_minimizers_all: usize = 0usize;
    #[cfg(debug_assertions)]
    let mut _total_query_owned_bytes_all: usize = 0usize;
    #[cfg(debug_assertions)]
    let mut max_query_owned_bytes: usize = 0usize;
    #[cfg(debug_assertions)]
    let mut max_mapping_result_bytes: usize = 0usize;
    let mut mapping_count: usize = 0usize;
    #[cfg(debug_assertions)]
    let mut mapping_detail_metrics: MappingMetrics = MappingMetrics::default();
    let mut emitted_pairs: usize = 0usize;
    let mut skipped_queries: usize = 0usize;
    let mut output: Box<dyn Write> = match &args.out_path {
        Some(path) => Box::new(BufWriter::new(fs::File::create(path)?)),
        None => Box::new(BufWriter::new(io::stdout())),
    };
    if args.emit_header {
        write_results_header(&mut *output)?;
    }
    let mut mapping_stats_output: Option<Box<dyn Write>> = match &args.mapping_stats_path {
        Some(path) => {
            let mut writer: Box<dyn Write> = Box::new(BufWriter::new(fs::File::create(path)?));
            write_mapping_stats_header(&mut *writer)?;
            Some(writer)
        }
        None => None,
    };

    if let SketchDatabase::Sharded { prefix, manifest } = &reference_database {
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
                                .map(|index| index.to_string())
                                .collect::<Vec<_>>()
                                .join(", "),
                            available
                                .iter()
                                .map(|index| index.to_string())
                                .collect::<Vec<_>>()
                                .join(", ")
                        ),
                    ));
                }
                manifest
                    .shards
                    .iter()
                    .cloned()
                    .filter(|shard| filter.contains(&shard.shard_index))
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

        let mut preloaded_queries: Vec<PreloadedQuery> = Vec::with_capacity(args.queries.len());
        for (query_index, query) in args.queries.iter().enumerate() {
            let query_path: String = query.label.clone();
            if progress_enabled {
                emit_progress(
                    "query",
                    &format!(
                        "event=start\tquery_done={query_index}\tquery_total={}\tpath={query_path}",
                        args.queries.len()
                    ),
                    total_start,
                );
            }

            let query_start: Instant = Instant::now();
            let mut reader_query: fasta::io::Reader<Box<dyn io::BufRead>> =
                open_fasta_reader(&query.open)?;
            let query_file: QueryFile = match QueryFile::collect(
                &mut reader_query,
                args.kmer_size,
                args.window_size,
                args.minmer_count,
                args.fragment_length,
                args.fragment_stride,
                args.min_fragment_length,
                args.split_n_run,
            ) {
                Ok(query_file) => query_file,
                Err(error) if args.queries.len() > 1 && is_no_usable_fragments_error(&error) => {
                    query_elapsed += query_start.elapsed();
                    skipped_queries += 1;
                    if progress_enabled {
                        eprintln!("WARNING\tskipping query with no usable fragments\t{query_path}");
                        emit_progress(
                            "query",
                            &format!(
                                "event=skipped\tquery_done={}\tquery_total={}\tpath={query_path}",
                                query_index + 1,
                                args.queries.len()
                            ),
                            total_start,
                        );
                    }
                    continue;
                }
                Err(error) => return Err(error),
            };
            query_elapsed += query_start.elapsed();

            let fragment_count: usize = query_file.fragments.len();
            let minimizer_count: usize = query_file.total_minimizers();
            let seed_minimizer_count: usize = query_file.total_seed_minimizers();
            total_query_fragments_all += fragment_count;
            total_query_minimizers_all += minimizer_count;
            total_query_seed_minimizers_all += seed_minimizer_count;

            #[cfg(debug_assertions)]
            {
                if performance_metrics_enabled {
                    let query_memory: QueryMemoryEstimate = query_file.memory_estimate();
                    let query_owned_bytes: usize = query_memory.fragment_struct_bytes
                        + query_memory.query_minimizer_vec_bytes
                        + query_memory.seed_minimizer_vec_bytes;
                    _total_query_owned_bytes_all += query_owned_bytes;
                    max_query_owned_bytes = max_query_owned_bytes.max(query_owned_bytes);
                }
            }

            preloaded_queries.push(PreloadedQuery {
                query_index,
                query_path,
                query_start,
                query_file,
                fragment_count,
                minimizer_count,
                seed_minimizer_count,
            });
        }

        if !preloaded_queries.is_empty() {
            let mut per_query_mapping_results: Vec<Vec<MappingResult>> =
                (0..preloaded_queries.len()).map(|_| Vec::new()).collect();
            let mut all_reference_files: Vec<ReferenceFile> =
                Vec::with_capacity(next_reference_file_offset);
            let mut all_reference_contig_names: Option<Vec<ReferenceContigName>> =
                mapping_stats_requested.then(|| Vec::with_capacity(next_reference_contig_offset));

            let (tx, rx) = mpsc::sync_channel::<io::Result<LoadedShard>>(1);
            let loader_active_shards = active_shards.clone();
            let loader_prefix = prefix.clone();
            let loader_tmp_dir = args.tmp_dir.clone();
            let loader_params = SketchParams {
                kmer_size: args.kmer_size,
                window_size: args.window_size,
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
                    .map(|sketch| {
                        sketch.prefetch_sequential();
                        LoadedShard {
                            shard_index: shard.shard_index,
                            shard_offset,
                            sketch,
                        }
                    });

                    if tx.send(load_result).is_err() {
                        break;
                    }
                }
            });

            let mut shards_done: usize = 0usize;
            while let Ok(load_result) = rx.recv() {
                let mut loaded: LoadedShard = load_result?;
                let frequency_threshold: usize = loaded
                    .sketch
                    .index
                    .frequency_threshold(args.freq_threshold_percent);
                let reference_file_offset: usize = reference_file_offsets[loaded.shard_offset];
                let reference_contig_offset: usize = reference_contig_offsets[loaded.shard_offset];
                let mut shard_mapping_count: usize = 0usize;

                for (query_slot, preloaded_query) in preloaded_queries.iter().enumerate() {
                    let mut raw_stats: RawQueryMappingStats = collect_query_mappings(
                        &loaded.sketch,
                        &preloaded_query.query_file,
                        args.kmer_size,
                        args.window_size,
                        args.mash_threshold,
                        args.mash_confidence,
                        args.threads,
                        frequency_threshold,
                        performance_metrics_enabled,
                    )?;
                    let query_shard_mapping_count: usize = raw_stats.mapping_results.len();
                    for mapping in &mut raw_stats.mapping_results {
                        mapping.reference_file_id += reference_file_offset;
                        mapping.reference_contig_id += reference_contig_offset;
                    }

                    mapping_elapsed += raw_stats.mapping_elapsed;
                    mapping_count += query_shard_mapping_count;
                    shard_mapping_count += query_shard_mapping_count;
                    #[cfg(debug_assertions)]
                    {
                        if performance_metrics_enabled {
                            mapping_detail_metrics.merge(raw_stats.metrics);
                            max_mapping_result_bytes =
                                max_mapping_result_bytes.max(raw_stats.max_mapping_result_bytes);
                        }
                    }
                    per_query_mapping_results[query_slot].append(&mut raw_stats.mapping_results);
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
            loader_handle.join().map_err(|_| {
                io::Error::new(io::ErrorKind::Other, "shard loader thread panicked")
            })?;

            for (query_slot, preloaded_query) in preloaded_queries.iter().enumerate() {
                let results: Vec<MappingResult> =
                    std::mem::take(&mut per_query_mapping_results[query_slot]);
                let stats = write_query_outputs(
                    &all_reference_files,
                    all_reference_contig_names.as_deref(),
                    results,
                    &preloaded_query.query_file,
                    &preloaded_query.query_path,
                    &mut *output,
                    mapping_stats_output.as_deref_mut(),
                    args.fragment_length,
                )?;
                summary_elapsed += stats.0;
                emitted_pairs += stats.1;

                let query_total_time = preloaded_query.query_start.elapsed();
                if progress_enabled {
                    emit_progress(
                        "query",
                        &format!(
                            "event=complete\tquery_done={}\tquery_total={}\tfragments={}\tminimizers={}\tseed_minimizers={}\telapsed_ms={}",
                            preloaded_query.query_index + 1,
                            args.queries.len(),
                            preloaded_query.fragment_count,
                            preloaded_query.minimizer_count,
                            preloaded_query.seed_minimizer_count,
                            query_total_time.as_millis()
                        ),
                        total_start,
                    );
                }
            }
        }
    }

    if let SketchDatabase::Single(reference_sketch) = &reference_database {
        for (query_index, query) in args.queries.iter().enumerate() {
            let query_path: &str = query.label.as_str();
            if progress_enabled {
                emit_progress(
                    "query",
                    &format!(
                        "event=start\tquery_done={query_index}\tquery_total={}\tpath={query_path}",
                        args.queries.len()
                    ),
                    total_start,
                );
            }

            let query_start: Instant = Instant::now();
            let mut reader_query: fasta::io::Reader<Box<dyn io::BufRead>> =
                open_fasta_reader(&query.open)?;
            let query_file: QueryFile = match QueryFile::collect(
                &mut reader_query,
                args.kmer_size,
                args.window_size,
                args.minmer_count,
                args.fragment_length,
                args.fragment_stride,
                args.min_fragment_length,
                args.split_n_run,
            ) {
                Ok(query_file) => query_file,
                Err(error) if args.queries.len() > 1 && is_no_usable_fragments_error(&error) => {
                    query_elapsed += query_start.elapsed();
                    skipped_queries += 1;
                    if progress_enabled {
                        eprintln!("WARNING\tskipping query with no usable fragments\t{query_path}");
                        emit_progress(
                            "query",
                            &format!(
                                "event=skipped\tquery_done={}\tquery_total={}\tpath={query_path}",
                                query_index + 1,
                                args.queries.len()
                            ),
                            total_start,
                        );
                    }
                    continue;
                }
                Err(error) => return Err(error),
            };
            query_elapsed += query_start.elapsed();

            let query_fragment_count: usize = query_file.fragments.len();
            let query_minimizer_count: usize = query_file.total_minimizers();
            let query_seed_minimizer_count: usize = query_file.total_seed_minimizers();
            total_query_fragments_all += query_fragment_count;
            total_query_minimizers_all += query_minimizer_count;
            total_query_seed_minimizers_all += query_seed_minimizer_count;

            #[cfg(debug_assertions)]
            {
                if performance_metrics_enabled {
                    let query_memory: QueryMemoryEstimate = query_file.memory_estimate();
                    let query_owned_bytes: usize = query_memory.fragment_struct_bytes
                        + query_memory.query_minimizer_vec_bytes
                        + query_memory.seed_minimizer_vec_bytes;
                    _total_query_owned_bytes_all += query_owned_bytes;
                    max_query_owned_bytes = max_query_owned_bytes.max(query_owned_bytes);
                }
            }

            let frequency_threshold: usize = reference_sketch
                .index
                .frequency_threshold(args.freq_threshold_percent);
            let stats: QueryMappingStats = map_query_against_reference_sketch(
                reference_sketch,
                &query_file,
                query_path,
                &args,
                frequency_threshold,
                performance_metrics_enabled,
                &mut *output,
                mapping_stats_output.as_deref_mut(),
            )?;
            mapping_elapsed += stats.mapping_elapsed;
            summary_elapsed += stats.summary_elapsed;
            mapping_count += stats.mapping_count;
            emitted_pairs += stats.emitted_pairs;
            #[cfg(debug_assertions)]
            {
                mapping_detail_metrics.merge(stats.metrics);
                max_mapping_result_bytes =
                    max_mapping_result_bytes.max(stats.max_mapping_result_bytes);
            }

            let query_total_time = query_start.elapsed();
            if progress_enabled {
                emit_progress(
                    "query",
                    &format!(
                        "event=complete\tquery_done={}\tquery_total={}\tfragments={query_fragment_count}\tminimizers={query_minimizer_count}\tseed_minimizers={query_seed_minimizer_count}\telapsed_ms={}",
                        query_index + 1,
                        args.queries.len(),
                        query_total_time.as_millis()
                    ),
                    total_start,
                );
            }
        }
    }

    output.flush()?;
    if let Some(mut stats_out) = mapping_stats_output {
        stats_out.flush()?;
    }

    if progress_enabled {
        emit_progress(
            "complete",
            &format!(
                "event=run\treference_mode={reference_mode}\treferences={}\tqueries={}\tskipped_queries={skipped_queries}\treference_contigs={}\tunique_minimizers={}\tquery_fragments={total_query_fragments_all}\tquery_minimizers={total_query_minimizers_all}\tseed_minimizers={total_query_seed_minimizers_all}\tmappings={mapping_count}\temitted_pairs={emitted_pairs}",
                reference_database.reference_count(),
                args.queries.len(),
                reference_database.contig_count(),
                reference_database.unique_minimizer_count()
            ),
            total_start,
        );
    }

    #[cfg(not(debug_assertions))]
    let _ = (
        reference_elapsed,
        query_elapsed,
        mapping_elapsed,
        summary_elapsed,
        performance_metrics_enabled,
    );

    #[cfg(debug_assertions)]
    {
        if performance_metrics_enabled {
            eprintln!(
                "METRICS\treference_mode={reference_mode}\tthreads={}\tfreq_threshold_percent={:.6}\tfreq_threshold={}\tminmer_count={}\tfragment_stride={}\tmin_fragment_length={}\tsplit_n_run={}\treferences={}\tqueries={}\tskipped_queries={skipped_queries}\treference_contigs={}\tunique_minimizers={}\tquery_fragments={total_query_fragments_all}\tmappings={mapping_count}\temitted_pairs={emitted_pairs}",
                args.threads,
                args.freq_threshold_percent,
                frequency_threshold_report,
                args.minmer_count
                    .map(|count| count.to_string())
                    .unwrap_or_else(|| "disabled".to_owned()),
                args.fragment_stride,
                args.min_fragment_length,
                args.split_n_run,
                reference_database.reference_count(),
                args.queries.len(),
                reference_database.contig_count(),
                reference_database.unique_minimizer_count(),
            );
            eprintln!(
                "METRICS\treference_ms={:.3}\tquery_ms={:.3}\tmapping_ms={:.3}\tsummary_ms={:.3}\ttotal_ms={:.3}\tpeak_rss_kb={}",
                reference_elapsed.as_secs_f64() * 1000.0,
                query_elapsed.as_secs_f64() * 1000.0,
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
                let upper_bound: String = SEED_HIT_HISTOGRAM_UPPER_BOUNDS
                    .get(bin_index)
                    .map(|upper_bound| upper_bound.to_string())
                    .unwrap_or_else(|| SEED_HIT_HISTOGRAM_OVERFLOW_LABEL.to_owned());
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
        }
    }

    Ok(())
}
