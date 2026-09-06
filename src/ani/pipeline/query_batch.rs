//! Query collection, scratch spooling, batching, and mapping accounting.

#[cfg(debug_assertions)]
use std::mem::size_of;
use std::{
    fs, io,
    io::{BufReader, BufWriter, Write},
    time::Instant,
};

use noodles::fasta;

use crate::ani::{
    cli::CliArgs,
    io_util::{open_fasta_reader, ScratchFile},
    mapping::{map_query_batch_to_reference_parallel, MappingExecutor},
    metrics::{MappingMetrics, MappingOutput},
    minimizer::is_no_usable_fragments_error,
    model::{
        query::{MappingResult, QueryFile},
        reference::ReferenceSketch,
    },
    runtime::emit_progress,
};

pub(super) struct RawQueryMappingStats {
    pub(crate) mapping_results: Vec<MappingResult>,
    pub(crate) metrics: MappingMetrics,
    #[cfg(debug_assertions)]
    pub(crate) max_mapping_result_bytes: usize,
}

pub(super) struct RawQueryMappingBatch {
    pub(crate) queries: Vec<RawQueryMappingStats>,
    pub(crate) mapping_elapsed: std::time::Duration,
}

pub(super) struct PreparedQuery {
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

pub(super) struct CollectedQuery {
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

pub(super) enum QueryPreparation {
    Ready(CollectedQuery),
    Skipped {
        collection_elapsed: std::time::Duration,
    },
}

#[derive(Default)]
pub(super) struct QueryAccounting {
    pub(crate) fragments: usize,
    pub(crate) minimizers: usize,
    pub(crate) seed_minimizers: usize,
    #[cfg(debug_assertions)]
    pub(crate) total_owned_bytes: usize,
    #[cfg(debug_assertions)]
    pub(crate) max_owned_bytes: usize,
}

impl QueryAccounting {
    pub(super) fn record(&mut self, query: &CollectedQuery, collect_memory_metrics: bool) {
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

pub(super) const MAX_QUERY_BATCH_SPOOL_BYTES: u64 = 256 * 1024 * 1024;
pub(super) const MAX_QUERY_BATCH_OWNED_BYTES: usize = 256 * 1024 * 1024;

pub(super) fn can_admit_query_owned_bytes(
    query_count: usize,
    current_owned_bytes: usize,
    candidate_owned_bytes: usize,
) -> bool {
    query_count == 0
        || current_owned_bytes.saturating_add(candidate_owned_bytes) <= MAX_QUERY_BATCH_OWNED_BYTES
}

impl PreparedQuery {
    pub(super) fn load(&self) -> io::Result<QueryFile> {
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

pub(super) fn prepare_query(
    args: &CliArgs,
    query_index: usize,
    completed_queries: &mut usize,
    total_start: Instant,
    progress_enabled: bool,
) -> io::Result<QueryPreparation> {
    let query = &args.queries[query_index];
    let query_path: String = query.output_label.clone();
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
    let mut reader_query: fasta::io::Reader<Box<dyn io::BufRead>> =
        open_fasta_reader(&query.input_path)?;
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

pub(super) fn spool_query_file(
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
pub(super) fn collect_query_mapping_batch(
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
