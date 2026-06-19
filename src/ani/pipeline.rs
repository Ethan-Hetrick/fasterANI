//! Top-level `run` orchestration and output writing.

#[cfg(debug_assertions)]
use std::mem::size_of;
use std::{
    collections::HashSet,
    fs, io,
    io::{BufWriter, Write},
    sync::atomic::{AtomicUsize, Ordering as AtomicOrdering},
    time::Instant,
};

use noodles::fasta;
use rayon::prelude::*;

use crate::ani::{
    check_memory_limit, emit_progress, fastani_compatible_fragment_mode, final_ani_computation,
    is_no_usable_fragments_error, legacy_sketch_path, manifest_path,
    map_query_to_reference_parallel, open_fasta_reader, parse_cli_args,
    performance_metrics_enabled, shard_entry_path, AniComputation, CliArgs, MappingOutput,
    MappingResult, MappingResultKey, QueryFile, QueryFragment, ReferenceContigName, ReferenceFile,
    ReferenceSketch, RuntimeOptions, ShardedBuildOptions, SketchDatabase, SketchParams,
};
#[cfg(debug_assertions)]
use crate::ani::{
    memory_gib, memory_mib, peak_rss_kb, MappingMetrics, MappingScratch, MinimizerKey,
    QueryMemoryEstimate, ReferenceMemoryEstimate, ReferenceMinimizer, SeedHit,
    SEED_HIT_HISTOGRAM_OVERFLOW_LABEL, SEED_HIT_HISTOGRAM_UPPER_BOUNDS, SKETCH_KEY_MODE,
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

struct ShardQueryResult {
    pub(crate) shard_index: usize,
    pub(crate) reference_files: Vec<ReferenceFile>,
    pub(crate) reference_contig_names: Option<Vec<ReferenceContigName>>,
    pub(crate) mapping_results: Vec<MappingResult>,
    pub(crate) load_elapsed: std::time::Duration,
    pub(crate) mapping_elapsed: std::time::Duration,
    pub(crate) mapping_count: usize,
    #[cfg(debug_assertions)]
    pub(crate) max_mapping_result_bytes: usize,
    #[cfg(debug_assertions)]
    pub(crate) metrics: MappingMetrics,
}

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
    // shared_/total_fragment_equivalents are fractional fragment counts
    // (aligned bases / fragment-length), matching FastANI's column layout, so
    // they are not necessarily whole numbers despite the `.2` formatting.
    writeln!(
        output,
        "query_file\treference_file\tani\tshared_fragment_equivalents\ttotal_fragment_equivalents"
    )
}

fn write_mapping_stats_header(output: &mut dyn Write) -> io::Result<()> {
    // Coordinates are 0-based half-open so the start/end columns can be used as BED intervals.
    writeln!(
        output,
        "query_file\treference_file\tquery_contig\treference_contig\tquery_fragment_id\tquery_start\tquery_end\treference_start\treference_end\tidentity\tquery_minimizer_count\treference_minimizer_count\tshared_minimizers\tunion_minimizers\tjaccard\tfragment_length\tis_reciprocal_best"
    )
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
        let reference_contig: String = if let Some(contig_name) =
            reference_contig_names.and_then(|contigs| contigs.get(mapping.reference_contig_id))
        {
            contig_name.name.clone()
        } else {
            mapping.reference_contig_id.to_string()
        };
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

        writeln!(
            output,
            "{query_path}\t{}\t{query_contig}\t{reference_contig}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{}\t{}\t{}\t{}\t{:.8}\t{}\t{}",
            reference_file.path,
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
    let summary_start: Instant = Instant::now();
    let mapping_results_for_stats: Option<Vec<MappingResult>> = mapping_stats_output
        .as_ref()
        .map(|_| mapping_results.clone());
    let ani_computation: AniComputation =
        final_ani_computation(mapping_results, reference_files.len(), fragment_length);
    let summary_elapsed: std::time::Duration = summary_start.elapsed();

    if let (Some(stats_results), Some(stats_output)) =
        (mapping_results_for_stats.as_ref(), mapping_stats_output)
    {
        write_mapping_stats(
            reference_files,
            reference_contig_names,
            stats_results,
            &ani_computation.reciprocal_best_keys,
            query_file,
            query_path,
            stats_output,
        )?;
    }

    let query_mapped_length: u64 = query_file.mapped_length();
    let mut emitted_pairs: usize = 0usize;

    for (reference_file, summary) in reference_files.iter().zip(ani_computation.summaries) {
        if summary.shared_fragments == 0 {
            continue;
        }

        let ani: f64 = summary.weighted_identity_sum / summary.shared_bases as f64;
        // Columns 4/5 are "fragment equivalents": aligned (and total) bases expressed
        // in units of fragment-length, mirroring FastANI's matched/total fragment
        // counts. They are ratios, so they can be fractional even though printed `.2`.
        let shared_fragment_equivalents: f64 = summary.shared_bases as f64 / fragment_length as f64;
        let total_fragment_equivalents: f64 = query_mapped_length as f64 / fragment_length as f64;
        writeln!(
            output,
            "{query_path}\t{}\t{ani:.3}\t{shared_fragment_equivalents:.2}\t{total_fragment_equivalents:.2}",
            reference_file.path
        )?;
        emitted_pairs += 1;
    }

    Ok((summary_elapsed, emitted_pairs))
}

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
        args.min_identity,
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
        max_memory_bytes: args.max_memory_bytes,
        worker_threads: args.threads,
    };
    check_memory_limit("startup", runtime_options)?;
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
                "event=algorithm\tkmer_size={}\twindow_size={}\tfragment_length={}\tfragment_stride={}\tmin_fragment_length={}\tmin_identity={:.6}\tmash_confidence={:.6}\tminmer_count={}\tfreq_threshold_percent={:.6}\tsplit_n_run={}",
                args.kmer_size,
                args.window_size,
                args.fragment_length,
                args.fragment_stride,
                args.min_fragment_length,
                args.min_identity,
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
            shard_size: args.shard_size,
            shard_minimizers: args.shard_minimizers,
            index_build_mode: args.index_build_mode,
            threads: args.threads,
        },
        mapping_stats_requested,
        runtime_options,
    )?;
    let reference_mode: &str =
        reference_database.mode_name(sketch_was_requested, existing_database_loaded);
    let mut reference_elapsed: std::time::Duration = reference_start.elapsed();
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
    let mut total_query_owned_bytes_all: usize = 0usize;
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

    for (query_index, query_path) in args.queries.iter().enumerate() {
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
        check_memory_limit("before query", runtime_options)?;

        let query_start: Instant = Instant::now();
        let mut reader_query: fasta::io::Reader<Box<dyn io::BufRead>> =
            open_fasta_reader(query_path)?;
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
                check_memory_limit("after skipped query", runtime_options)?;
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
                total_query_owned_bytes_all += query_owned_bytes;
                max_query_owned_bytes = max_query_owned_bytes.max(query_owned_bytes);
            }
        }

        match &reference_database {
            SketchDatabase::Single(reference_sketch) => {
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
            }
            SketchDatabase::Sharded { prefix, manifest } => {
                let mut query_reference_files: Vec<ReferenceFile> =
                    Vec::with_capacity(manifest.total_references);
                let mut query_reference_contig_names: Option<Vec<ReferenceContigName>> =
                    mapping_stats_requested
                        .then(|| Vec::with_capacity(manifest.total_reference_contigs));
                let mut query_mapping_results: Vec<MappingResult> = Vec::new();
                let shard_query_parallelism: usize =
                    args.threads.max(1).min(manifest.shards.len().max(1));
                let mapping_threads_per_shard: usize =
                    args.threads.max(1).div_ceil(shard_query_parallelism.max(1));
                let mut reference_contig_offsets: Vec<usize> =
                    Vec::with_capacity(manifest.shards.len());
                let mut next_reference_contig_offset: usize = 0usize;
                for shard in &manifest.shards {
                    reference_contig_offsets.push(next_reference_contig_offset);
                    next_reference_contig_offset += shard.reference_contigs;
                }
                let completed_shard_queries: AtomicUsize = AtomicUsize::new(0);
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(shard_query_parallelism)
                    .build()
                    .map_err(|err| {
                        io::Error::new(
                            io::ErrorKind::InvalidInput,
                            format!("failed to initialize sharded query thread pool: {err}"),
                        )
                    })?;
                let mut shard_results: Vec<ShardQueryResult> = pool.install(|| {
                    manifest
                        .shards
                        .par_iter()
                        .enumerate()
                        .map(|(shard_offset, shard)| {
                            if progress_enabled {
                                emit_progress(
                                    "shard_load",
                                    &format!(
                                        "event=start\tquery_done={query_index}\tshard={}\tshards_total={}\tfilename={}\tshard_query_parallelism={shard_query_parallelism}\tmapping_threads_per_shard={mapping_threads_per_shard}",
                                        shard.shard_index,
                                        manifest.shards.len(),
                                        shard.filename
                                    ),
                                    total_start,
                                );
                            }
                            let shard_load_start: Instant = Instant::now();
                            let shard_sketch: ReferenceSketch = ReferenceSketch::load(
                                &shard_entry_path(prefix, shard),
                                SketchParams {
                                    kmer_size: args.kmer_size,
                                    window_size: args.window_size,
                                    fragment_length: args.fragment_length,
                                    min_fragment_length: args.min_fragment_length,
                                    split_n_run: args.split_n_run,
                                },
                                mapping_stats_requested,
                                args.tmp_dir.as_deref(),
                                runtime_options,
                            )?;
                            let load_elapsed: std::time::Duration = shard_load_start.elapsed();
                            let frequency_threshold: usize = shard_sketch
                                .index
                                .frequency_threshold(args.freq_threshold_percent);
                            let mut raw_stats: RawQueryMappingStats = collect_query_mappings(
                                &shard_sketch,
                                &query_file,
                                args.kmer_size,
                                args.window_size,
                                args.min_identity,
                                args.mash_confidence,
                                mapping_threads_per_shard,
                                frequency_threshold,
                                performance_metrics_enabled,
                            )?;
                            let shard_mapping_count: usize = raw_stats.mapping_results.len();
                            let reference_file_offset: usize = shard.first_reference;
                            let reference_contig_offset: usize =
                                reference_contig_offsets[shard_offset];

                            for mapping in &mut raw_stats.mapping_results {
                                mapping.reference_file_id += reference_file_offset;
                                mapping.reference_contig_id += reference_contig_offset;
                            }
                            let reference_files: Vec<ReferenceFile> =
                                shard_sketch.files.to_vec();
                            let reference_contig_names: Option<Vec<ReferenceContigName>> =
                                shard_sketch.contig_names.clone();
                            let shards_done: usize =
                                completed_shard_queries.fetch_add(1, AtomicOrdering::Relaxed) + 1;

                            if progress_enabled {
                                emit_progress(
                                    "shard_load",
                                    &format!(
                                        "event=complete\tquery_done={query_index}\tshard={}\tshards_done={shards_done}\tshards_total={}\tmappings={shard_mapping_count}",
                                        shard.shard_index,
                                        manifest.shards.len()
                                    ),
                                    total_start,
                                );
                            }
                            drop(shard_sketch);
                            check_memory_limit("after shard query", runtime_options)?;

                            Ok(ShardQueryResult {
                                shard_index: shard.shard_index,
                                reference_files,
                                reference_contig_names,
                                mapping_results: raw_stats.mapping_results,
                                load_elapsed,
                                mapping_elapsed: raw_stats.mapping_elapsed,
                                mapping_count: shard_mapping_count,
                                #[cfg(debug_assertions)]
                                max_mapping_result_bytes: raw_stats.max_mapping_result_bytes,
                                #[cfg(debug_assertions)]
                                metrics: raw_stats.metrics,
                            })
                        })
                        .collect::<io::Result<Vec<_>>>()
                })?;
                shard_results.sort_by_key(|result| result.shard_index);

                for mut shard_result in shard_results {
                    reference_elapsed += shard_result.load_elapsed;
                    mapping_elapsed += shard_result.mapping_elapsed;
                    mapping_count += shard_result.mapping_count;
                    query_reference_files.append(&mut shard_result.reference_files);
                    if let Some(all_contig_names) = query_reference_contig_names.as_mut() {
                        let mut shard_contig_names: Vec<ReferenceContigName> =
                            shard_result.reference_contig_names.ok_or_else(|| {
                                io::Error::new(
                                    io::ErrorKind::InvalidData,
                                    "reference contig names were not loaded for mapping stats",
                                )
                            })?;
                        all_contig_names.append(&mut shard_contig_names);
                    }
                    query_mapping_results.append(&mut shard_result.mapping_results);
                    #[cfg(debug_assertions)]
                    {
                        mapping_detail_metrics.merge(shard_result.metrics);
                        max_mapping_result_bytes =
                            max_mapping_result_bytes.max(shard_result.max_mapping_result_bytes);
                    }
                }

                let (query_summary_elapsed, query_emitted_pairs) = write_query_outputs(
                    &query_reference_files,
                    query_reference_contig_names.as_deref(),
                    query_mapping_results,
                    &query_file,
                    query_path,
                    &mut *output,
                    mapping_stats_output.as_deref_mut(),
                    args.fragment_length,
                )?;
                summary_elapsed += query_summary_elapsed;
                emitted_pairs += query_emitted_pairs;
            }
        }

        if progress_enabled {
            emit_progress(
                "query",
                &format!(
                    "event=complete\tquery_done={}\tquery_total={}\tpath={query_path}\tfragments={query_fragment_count}\tquery_minimizers={query_minimizer_count}\tseed_minimizers={query_seed_minimizer_count}\tmappings={}",
                    query_index + 1,
                    args.queries.len(),
                    mapping_count
                ),
                total_start,
            );
        }
        check_memory_limit("after query", runtime_options)?;
    }
    output.flush()?;

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
    if performance_metrics_enabled {
        eprintln!(
            "METRICS\treference_mode={reference_mode}\tthreads={}\tmax_memory_gb={}\tfreq_threshold_percent={:.6}\tfreq_threshold={}\tminmer_count={}\tfragment_stride={}\tmin_fragment_length={}\tsplit_n_run={}\treferences={}\tqueries={}\tskipped_queries={skipped_queries}\treference_contigs={}\tunique_minimizers={}\tquery_fragments={total_query_fragments_all}\tmappings={mapping_count}\temitted_pairs={emitted_pairs}",
            args.threads,
            args.max_memory_bytes
                .map(|bytes| format!("{:.3}", memory_gib(bytes)))
                .unwrap_or_else(|| "disabled".to_owned()),
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
            mapping_detail_metrics
                .candidate_discovery_elapsed
                .as_secs_f64()
                * 1000.0,
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
            eprintln!(
                "MEMORY_ESTIMATE\tscope=mmap_arrays\tslot_key_bytes={}\tslot_key_mib={:.3}\thit_offset_bytes={}\thit_offset_mib={:.3}\thit_count_bytes={}\thit_count_mib={:.3}\thit_payload_bytes={}\thit_payload_mib={:.3}\tcontig_record_bytes={}\tcontig_record_mib={:.3}\treference_minimizer_bytes={}\treference_minimizer_mib={:.3}",
                reference_memory.mmap_slot_key_bytes,
                memory_mib(reference_memory.mmap_slot_key_bytes),
                reference_memory.mmap_hit_offset_bytes,
                memory_mib(reference_memory.mmap_hit_offset_bytes),
                reference_memory.mmap_hit_count_bytes,
                memory_mib(reference_memory.mmap_hit_count_bytes),
                reference_memory.mmap_hit_payload_bytes,
                memory_mib(reference_memory.mmap_hit_payload_bytes),
                reference_memory.mmap_contig_record_bytes,
                memory_mib(reference_memory.mmap_contig_record_bytes),
                reference_memory.mmap_reference_minimizer_bytes,
                memory_mib(reference_memory.mmap_reference_minimizer_bytes),
            );
        } else if let SketchDatabase::Sharded { manifest, .. } = &reference_database {
            eprintln!(
                "MEMORY_ESTIMATE\tscope=reference_sharded\tshards={}\treferences={}\treference_contigs={}\treference_minimizers={}\tunique_minimizers_shard_sum={}\tmapped_reference_length={}",
                manifest.shards.len(),
                manifest.total_references,
                manifest.total_reference_contigs,
                manifest.total_reference_minimizers,
                manifest.total_shard_unique_minimizers,
                manifest.total_mapped_reference_length,
            );
        }
        eprintln!(
            "MEMORY_ESTIMATE\tscope=query_total\tquery_fragments={}\tquery_minimizers={}\tquery_seed_minimizers={}\tquery_owned_bytes_sum={}\tquery_owned_mib_sum={:.3}\tmax_query_owned_bytes={}\tmax_query_owned_mib={:.3}",
            total_query_fragments_all,
            total_query_minimizers_all,
            total_query_seed_minimizers_all,
            total_query_owned_bytes_all,
            memory_mib(total_query_owned_bytes_all),
            max_query_owned_bytes,
            memory_mib(max_query_owned_bytes),
        );
        eprintln!(
            "MEMORY_ESTIMATE\tscope=mapping\tmapping_results={mapping_count}\tmax_mapping_result_bytes={}\tmax_mapping_result_mib={:.3}\tthreads={}\tmapping_scratch_struct_only_bytes={}\tmapping_scratch_struct_only_mib={:.3}",
            max_mapping_result_bytes,
            memory_mib(max_mapping_result_bytes),
            args.threads,
            args.threads * size_of::<MappingScratch>(),
            memory_mib(args.threads * size_of::<MappingScratch>()),
        );

        let reference_minimizer_count_for_projection: usize = match &reference_database {
            SketchDatabase::Single(_) => reference_memory
                .as_ref()
                .map(|memory| memory.reference_minimizers)
                .unwrap_or(0),
            SketchDatabase::Sharded { manifest, .. } => manifest.total_reference_minimizers,
        };
        let unique_index_keys_for_projection: usize = match &reference_database {
            SketchDatabase::Single(_) => reference_memory
                .as_ref()
                .map(|memory| memory.unique_index_keys)
                .unwrap_or(0),
            SketchDatabase::Sharded { manifest, .. } => manifest.total_shard_unique_minimizers,
        };
        let reference_minimizer_u32_savings: usize = reference_minimizer_count_for_projection
            .saturating_mul(size_of::<ReferenceMinimizer>().saturating_sub(size_of::<SeedHit>()));
        let index_key_u32_savings: usize = unique_index_keys_for_projection
            .saturating_mul(size_of::<MinimizerKey>().saturating_sub(size_of::<u32>()));
        let query_key_u32_savings: usize = total_query_minimizers_all
            .saturating_add(total_query_seed_minimizers_all)
            .saturating_mul(size_of::<MinimizerKey>().saturating_sub(size_of::<u32>()));
        let total_u32_projected_savings: usize = reference_minimizer_u32_savings
            .saturating_add(index_key_u32_savings)
            .saturating_add(query_key_u32_savings);
        let peak_rss_bytes: usize = match peak_rss_kb() {
            peak if peak > 0 => peak as usize * 1024,
            _ => 0,
        };
        let projected_savings_peak_percent: f64 = if peak_rss_bytes == 0 {
            f64::NAN
        } else {
            100.0 * total_u32_projected_savings as f64 / peak_rss_bytes as f64
        };
        eprintln!(
            "MEMORY_ESTIMATE\tscope=u32_key_projection\tkey_mode={SKETCH_KEY_MODE}\treference_minimizer_savings_bytes={}\tindex_key_savings_bytes={}\tquery_key_savings_bytes={}\ttotal_projected_savings_bytes={}\ttotal_projected_savings_mib={:.3}\tprojected_savings_vs_peak_rss_percent={:.3}",
            reference_minimizer_u32_savings,
            index_key_u32_savings,
            query_key_u32_savings,
            total_u32_projected_savings,
            memory_mib(total_u32_projected_savings),
            projected_savings_peak_percent,
        );
    }

    Ok(())
}
