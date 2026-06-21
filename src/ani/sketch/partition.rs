//! Index build modes and shard/partition planning.

use std::{
    fs, io,
    io::{BufWriter, Write},
    mem::size_of,
    path::Path,
    sync::atomic::{AtomicUsize, Ordering as AtomicOrdering},
    time::Instant,
};

use noodles::fasta;
use rayon::prelude::*;

use crate::ani::{
    check_memory_limit, emit_progress, expected_minimizer_window_count, open_fasta_reader,
    slice_as_bytes, split_sequence_ranges, validate_shard_minimizers, validate_shard_size,
    FastaInput, MinimizerKey, RuntimeOptions, ScratchFile, SeedHit, ShardManifest, ShardPlan,
    DEFAULT_PARTITION_TARGET_BYTES, ESTIMATED_PARTITIONED_SHARD_BYTES_PER_MINIMIZER,
    MAX_PARTITION_COUNT, MIN_PARTITION_COUNT, PARTITIONED_INDEX_MINIMIZER_THRESHOLD,
    REFERENCE_PROGRESS_INTERVAL, SKETCH_DATABASE_SCHEMA_VERSION, SKETCH_KEY_MODE, SKETCH_VERSION,
};

/// User-selectable strategy for building the reference sketch index.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum IndexBuildMode {
    Auto,
    Hash,
    Partitioned,
}

impl IndexBuildMode {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Auto => "auto",
            Self::Hash => "hash",
            Self::Partitioned => "partitioned",
        }
    }

    pub(crate) fn parse(value: &str) -> io::Result<Self> {
        match value {
            "auto" => Ok(Self::Auto),
            "hash" => Ok(Self::Hash),
            "partitioned" => Ok(Self::Partitioned),
            _ => Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("--index-build-mode must be one of auto, hash, partitioned; got {value:?}"),
            )),
        }
    }
}

/// One ungrouped minimizer-index hit written to a temporary partition.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(crate) struct PartitionHitRecord {
    pub(crate) key: MinimizerKey,
    pub(crate) hit: SeedHit,
}

/// One grouped minimizer key and the hit range assigned to it in the final payload file.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(crate) struct GroupedKeyRecord {
    pub(crate) key: MinimizerKey,
    pub(crate) hit_offset: u64,
    pub(crate) hit_count: u32,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct PartitionBuildPlan {
    pub(crate) partition_count: usize,
    pub(crate) target_partition_bytes: usize,
    pub(crate) estimated_record_bytes: usize,
}

pub(crate) struct PartitionGroupResult {
    pub(crate) partition_index: usize,
    pub(crate) grouped_key_scratch: ScratchFile,
    pub(crate) hit_payload_scratch: ScratchFile,
    pub(crate) keys: Vec<MinimizerKey>,
    pub(crate) key_count: usize,
    pub(crate) hit_count: usize,
}

pub(crate) fn estimate_partitioned_shard_memory_bytes(estimated_minimizers: usize) -> usize {
    estimated_minimizers
        .saturating_mul(ESTIMATED_PARTITIONED_SHARD_BYTES_PER_MINIMIZER)
        .max(DEFAULT_PARTITION_TARGET_BYTES)
}

pub(crate) fn database_build_parallelism(
    threads: usize,
    shard_plans: &[ShardPlan],
    max_memory_bytes: Option<u64>,
) -> usize {
    let requested_threads: usize = threads.max(1);
    let shard_count: usize = shard_plans.len().max(1);
    let thread_limited_jobs: usize = requested_threads.min(shard_count);

    let Some(max_memory_bytes) = max_memory_bytes else {
        return thread_limited_jobs;
    };

    let Ok(max_memory_bytes) = usize::try_from(max_memory_bytes) else {
        return thread_limited_jobs;
    };

    let largest_shard_memory_bytes: usize = shard_plans
        .iter()
        .map(|plan| estimate_partitioned_shard_memory_bytes(plan.estimated_minimizers))
        .max()
        .unwrap_or(DEFAULT_PARTITION_TARGET_BYTES);
    let memory_limited_jobs: usize = (max_memory_bytes / largest_shard_memory_bytes).max(1);

    thread_limited_jobs.min(memory_limited_jobs)
}

fn ceil_div_usize(numerator: usize, denominator: usize) -> usize {
    if denominator == 0 {
        return numerator;
    }

    numerator.saturating_add(denominator.saturating_sub(1)) / denominator
}

pub(crate) fn partition_build_plan(
    estimated_minimizers: usize,
    max_memory_bytes: Option<u64>,
) -> PartitionBuildPlan {
    let estimated_record_bytes: usize =
        estimated_minimizers.saturating_mul(size_of::<PartitionHitRecord>());
    let memory_target: usize = max_memory_bytes
        .and_then(|bytes| usize::try_from(bytes / 8).ok())
        .unwrap_or(DEFAULT_PARTITION_TARGET_BYTES);
    let target_partition_bytes: usize = DEFAULT_PARTITION_TARGET_BYTES
        .min(memory_target)
        .max(size_of::<PartitionHitRecord>());
    let raw_partition_count: usize =
        ceil_div_usize(estimated_record_bytes.max(1), target_partition_bytes);
    let next_power: usize = raw_partition_count
        .checked_next_power_of_two()
        .unwrap_or(MAX_PARTITION_COUNT);
    let partition_count: usize = next_power.clamp(MIN_PARTITION_COUNT, MAX_PARTITION_COUNT);

    PartitionBuildPlan {
        partition_count,
        target_partition_bytes,
        estimated_record_bytes,
    }
}

pub(crate) fn effective_index_build_mode(
    requested_mode: IndexBuildMode,
    estimated_minimizers: usize,
) -> IndexBuildMode {
    match requested_mode {
        IndexBuildMode::Auto if estimated_minimizers >= PARTITIONED_INDEX_MINIMIZER_THRESHOLD => {
            IndexBuildMode::Partitioned
        }
        IndexBuildMode::Auto => IndexBuildMode::Hash,
        mode => mode,
    }
}

pub(crate) fn partition_id_for_key(key: MinimizerKey, partition_count: usize) -> usize {
    debug_assert!(partition_count.is_power_of_two());
    let partition_bits: u32 = partition_count.trailing_zeros();
    if partition_bits == 0 {
        return 0;
    }

    (key >> (MinimizerKey::BITS - partition_bits)) as usize
}

pub(crate) fn estimate_reference_minimizer_windows(
    reference: &FastaInput,
    kmer_size: usize,
    window_size: usize,
    split_n_run: usize,
) -> io::Result<usize> {
    let mut reader: fasta::io::Reader<Box<dyn io::BufRead>> = open_fasta_reader(&reference.open)?;
    let mut minimizer_window_count: usize = 0usize;

    for result in reader.records() {
        let record: fasta::Record = result.map_err(|err| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "failed to read FASTA record from reference {}: {err}",
                    reference.label
                ),
            )
        })?;
        let sequence: &fasta::record::Sequence = record.sequence();
        let sequence_bytes: &[u8] = sequence.as_ref();

        for segment_range in split_sequence_ranges(sequence_bytes, split_n_run) {
            minimizer_window_count +=
                expected_minimizer_window_count(segment_range.len(), kmer_size, window_size);
        }
    }

    Ok(minimizer_window_count)
}

pub(crate) fn estimate_selected_minimizers_from_windows(
    window_count: usize,
    window_size: usize,
) -> usize {
    let denominator: usize = window_size.saturating_add(1);
    if denominator == 0 {
        return window_count;
    }

    window_count
        .saturating_mul(2)
        .saturating_add(denominator - 1)
        / denominator
}

pub(crate) fn plan_shards_from_minimizer_counts(
    minimizer_counts: &[usize],
    shard_size: usize,
    shard_minimizers: usize,
) -> io::Result<Vec<ShardPlan>> {
    validate_shard_size(shard_size)?;
    validate_shard_minimizers(shard_minimizers)?;

    let mut plans: Vec<ShardPlan> = Vec::new();
    let mut shard_first_reference: usize = 0usize;
    let mut shard_reference_count: usize = 0usize;
    let mut shard_estimated_minimizers: usize = 0usize;

    for (reference_index, reference_minimizers) in minimizer_counts.iter().copied().enumerate() {
        let would_exceed_minimizers: bool = shard_reference_count > 0
            && shard_estimated_minimizers.saturating_add(reference_minimizers) > shard_minimizers;
        let would_exceed_reference_count: bool = shard_reference_count >= shard_size;

        if would_exceed_minimizers || would_exceed_reference_count {
            plans.push(ShardPlan {
                first_reference: shard_first_reference,
                reference_count: shard_reference_count,
                estimated_minimizers: shard_estimated_minimizers,
            });
            shard_first_reference = reference_index;
            shard_reference_count = 0;
            shard_estimated_minimizers = 0;
        }

        shard_reference_count += 1;
        shard_estimated_minimizers =
            shard_estimated_minimizers.saturating_add(reference_minimizers);
    }

    if shard_reference_count > 0 {
        plans.push(ShardPlan {
            first_reference: shard_first_reference,
            reference_count: shard_reference_count,
            estimated_minimizers: shard_estimated_minimizers,
        });
    }

    Ok(plans)
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn plan_shards_by_minimizers(
    references: &[FastaInput],
    kmer_size: usize,
    window_size: usize,
    split_n_run: usize,
    shard_size: usize,
    shard_minimizers: usize,
    threads: usize,
    runtime_options: RuntimeOptions,
) -> io::Result<Vec<ShardPlan>> {
    validate_shard_size(shard_size)?;
    validate_shard_minimizers(shard_minimizers)?;

    let plan_start: Instant = Instant::now();
    let planner_threads: usize = threads.max(1).min(references.len().max(1));

    if runtime_options.progress_enabled {
        emit_progress(
            "shard_plan",
            &format!(
                "event=start\testimator=window-count\treferences={}\tshard_size={shard_size}\tshard_minimizers={shard_minimizers}\tthreads={threads}\tplanner_parallelism={planner_threads}",
                references.len(),
            ),
            plan_start,
        );
    }

    let completed_references: AtomicUsize = AtomicUsize::new(0);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(planner_threads)
        .build()
        .map_err(|err| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("failed to initialize shard planning thread pool: {err}"),
            )
        })?;
    let minimizer_counts: Vec<usize> = pool.install(|| {
        references
            .par_iter()
            .map(|reference| {
                let reference_minimizer_windows: usize = estimate_reference_minimizer_windows(
                    reference,
                    kmer_size,
                    window_size,
                    split_n_run,
                )?;
                let reference_minimizers: usize = estimate_selected_minimizers_from_windows(
                    reference_minimizer_windows,
                    window_size,
                );

                let references_done: usize =
                    completed_references.fetch_add(1, AtomicOrdering::Relaxed) + 1;
                if runtime_options.progress_enabled
                    && (references_done.is_multiple_of(REFERENCE_PROGRESS_INTERVAL)
                        || references_done == references.len())
                {
                    emit_progress(
                        "shard_plan",
                        &format!(
                            "event=references\treferences_done={references_done}\treferences_total={}\tplanner_parallelism={planner_threads}",
                            references.len(),
                        ),
                        plan_start,
                    );
                }
                check_memory_limit("during shard planning", runtime_options)?;

                Ok(reference_minimizers)
            })
            .collect::<io::Result<Vec<usize>>>()
    })?;

    let plans: Vec<ShardPlan> =
        plan_shards_from_minimizer_counts(&minimizer_counts, shard_size, shard_minimizers)?;

    if runtime_options.progress_enabled {
        emit_progress(
            "shard_plan",
            &format!("event=complete\tshards={}", plans.len()),
            plan_start,
        );
    }

    Ok(plans)
}

pub(crate) fn shard_manifest_compatibility_error(
    manifest: &ShardManifest,
    kmer_size: usize,
    window_size: usize,
    fragment_length: u32,
    min_fragment_length: u32,
    split_n_run: usize,
) -> Option<String> {
    if manifest.dust_enabled {
        return Some(
            "reference sketch database is incompatible: it was built with the removed --dust filter"
                .to_owned(),
        );
    }

    if manifest.sketch_format_version != SKETCH_VERSION
        || manifest.database_schema_version != SKETCH_DATABASE_SCHEMA_VERSION
        || manifest.k != kmer_size
        || manifest.w != window_size
        || manifest.key_mode != SKETCH_KEY_MODE
        || manifest.fragment_length != fragment_length
        || manifest.min_fragment_length != min_fragment_length
        || manifest.split_n_run != split_n_run
    {
        return Some(format!(
            "reference sketch database is incompatible: sketch_format_version={} database_schema_version={} k={} w={} key_mode={} fragment_length={} min_fragment_length={} split_n_run={}",
            manifest.sketch_format_version,
            manifest.database_schema_version,
            manifest.k,
            manifest.w,
            manifest.key_mode,
            manifest.fragment_length,
            manifest.min_fragment_length,
            manifest.split_n_run
        ));
    }

    None
}

pub(crate) struct PartitionWriters {
    pub(crate) files: Vec<ScratchFile>,
    pub(crate) writers: Vec<BufWriter<fs::File>>,
    pub(crate) buffers: Vec<Vec<PartitionHitRecord>>,
    pub(crate) buffer_record_limit: usize,
}

impl PartitionWriters {
    pub(crate) fn new(
        partition_count: usize,
        tmp_dir: Option<&Path>,
        buffer_record_limit: usize,
    ) -> io::Result<Self> {
        if !partition_count.is_power_of_two() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "partition count must be a power of two",
            ));
        }

        let mut files: Vec<ScratchFile> = Vec::with_capacity(partition_count);
        let mut writers: Vec<BufWriter<fs::File>> = Vec::with_capacity(partition_count);
        for partition_index in 0..partition_count {
            let (scratch, file): (ScratchFile, fs::File) =
                ScratchFile::create(tmp_dir, &format!("index-partition-{partition_index}"))?;
            files.push(scratch);
            writers.push(BufWriter::new(file));
        }

        Ok(Self {
            files,
            writers,
            buffers: vec![Vec::new(); partition_count],
            buffer_record_limit: buffer_record_limit.max(1),
        })
    }

    pub(crate) fn partition_count(&self) -> usize {
        self.files.len()
    }

    pub(crate) fn path(&self, partition_index: usize) -> &Path {
        &self.files[partition_index].path
    }

    pub(crate) fn push(&mut self, record: PartitionHitRecord) -> io::Result<()> {
        let partition_index: usize = partition_id_for_key(record.key, self.partition_count());
        let buffer: &mut Vec<PartitionHitRecord> = &mut self.buffers[partition_index];
        buffer.push(record);

        if buffer.len() >= self.buffer_record_limit {
            self.flush_partition(partition_index)?;
        }

        Ok(())
    }

    pub(crate) fn flush_partition(&mut self, partition_index: usize) -> io::Result<()> {
        let buffer: &mut Vec<PartitionHitRecord> = &mut self.buffers[partition_index];
        if buffer.is_empty() {
            return Ok(());
        }

        self.writers[partition_index].write_all(slice_as_bytes(buffer))?;
        buffer.clear();

        Ok(())
    }

    pub(crate) fn flush_all(&mut self) -> io::Result<()> {
        for partition_index in 0..self.partition_count() {
            self.flush_partition(partition_index)?;
            self.writers[partition_index].flush()?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::ani::{
        database_build_parallelism, partition_build_plan, partition_id_for_key,
        plan_shards_from_minimizer_counts, PartitionBuildPlan, ShardPlan,
        DEFAULT_PARTITION_TARGET_BYTES, MAX_PARTITION_COUNT, MIN_PARTITION_COUNT,
    };
    use std::io;

    #[test]
    fn shard_planner_splits_by_minimizer_target() -> io::Result<()> {
        let plans: Vec<ShardPlan> =
            plan_shards_from_minimizer_counts(&[200, 250, 100, 400], 10, 500)?;

        assert_eq!(
            plans,
            vec![
                ShardPlan {
                    first_reference: 0,
                    reference_count: 2,
                    estimated_minimizers: 450,
                },
                ShardPlan {
                    first_reference: 2,
                    reference_count: 2,
                    estimated_minimizers: 500,
                },
            ]
        );

        Ok(())
    }

    #[test]
    fn shard_planner_respects_reference_count_cap() -> io::Result<()> {
        let plans: Vec<ShardPlan> =
            plan_shards_from_minimizer_counts(&[10, 10, 10, 10, 10], 2, 500)?;

        assert_eq!(
            plans,
            vec![
                ShardPlan {
                    first_reference: 0,
                    reference_count: 2,
                    estimated_minimizers: 20,
                },
                ShardPlan {
                    first_reference: 2,
                    reference_count: 2,
                    estimated_minimizers: 20,
                },
                ShardPlan {
                    first_reference: 4,
                    reference_count: 1,
                    estimated_minimizers: 10,
                },
            ]
        );

        Ok(())
    }

    #[test]
    fn automatic_partition_count_uses_power_of_two_bounds() {
        let small_plan: PartitionBuildPlan = partition_build_plan(10, None);
        assert_eq!(small_plan.partition_count, MIN_PARTITION_COUNT);
        assert!(small_plan.partition_count.is_power_of_two());

        let large_plan: PartitionBuildPlan = partition_build_plan(2_000_000_000, None);
        assert!(large_plan.partition_count > MIN_PARTITION_COUNT);
        assert!(large_plan.partition_count.is_power_of_two());
        assert!(large_plan.partition_count <= MAX_PARTITION_COUNT);

        let constrained_plan: PartitionBuildPlan =
            partition_build_plan(100_000_000, Some(1024 * 1024 * 1024));
        assert!(constrained_plan.target_partition_bytes < DEFAULT_PARTITION_TARGET_BYTES);
        assert!(constrained_plan.partition_count.is_power_of_two());
    }

    #[test]
    fn database_build_parallelism_uses_threads_when_memory_allows() {
        let shard_plans: Vec<ShardPlan> = vec![
            ShardPlan {
                first_reference: 0,
                reference_count: 1,
                estimated_minimizers: 300_000_000,
            };
            16
        ];
        let one_hundred_gib: u64 = 100 * 1024 * 1024 * 1024;

        assert_eq!(
            database_build_parallelism(12, &shard_plans, Some(one_hundred_gib)),
            12
        );
        assert_eq!(database_build_parallelism(4, &shard_plans, None), 4);
    }

    #[test]
    fn partition_id_uses_high_key_bits() {
        assert_eq!(partition_id_for_key(0x0000_0000, 16), 0);
        assert_eq!(partition_id_for_key(0x1000_0000, 16), 1);
        assert_eq!(partition_id_for_key(0xF000_0000, 16), 15);
    }
}
