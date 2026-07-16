//! Index build modes and shard/partition planning.

use std::{
    fs, io,
    io::{BufWriter, Write},
    mem::{offset_of, size_of},
    path::Path,
    str::FromStr,
    sync::atomic::{AtomicUsize, Ordering as AtomicOrdering},
    time::Instant,
};

use noodles::fasta;
use rayon::prelude::*;

use crate::ani::{
    emit_progress, expected_minimizer_window_count, open_fasta_reader, slice_as_bytes,
    split_sequence_ranges, validate_max_shard_minimizers, FastaInput, MinimizerKey, RuntimeOptions,
    ScratchFile, SeedHit, ShardManifest, ShardPlan, DEFAULT_PARTITION_TARGET_BYTES,
    ESTIMATED_PARTITIONED_SHARD_BYTES_PER_MINIMIZER, MAX_PARTITION_COUNT, MIN_PARTITION_COUNT,
    PARTITIONED_INDEX_MINIMIZER_THRESHOLD, REFERENCE_PROGRESS_INTERVAL,
    SKETCH_DATABASE_SCHEMA_VERSION, SKETCH_KEY_MODE, SKETCH_VERSION,
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
}

impl FromStr for IndexBuildMode {
    type Err = io::Error;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
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

// The AVX2 gather below uses these byte offsets directly. Keep the assumptions
// compile-time checked instead of relying on an undocumented repr(C) layout.
const _: [(); 12] = [(); size_of::<PartitionHitRecord>()];
const _: [(); 0] = [(); offset_of!(PartitionHitRecord, key)];
const _: [(); 4] = [(); offset_of!(PartitionHitRecord, hit)];
const _: [(); 8] = [(); size_of::<SeedHit>()];
const _: [(); 0] = [(); offset_of!(SeedHit, reference_contig_id)];
const _: [(); 4] = [(); offset_of!(SeedHit, position)];

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

pub(crate) fn database_build_parallelism(threads: usize, shard_plans: &[ShardPlan]) -> usize {
    let requested_threads: usize = threads.max(1);
    let shard_count: usize = shard_plans.len().max(1);
    let thread_limit: usize = requested_threads.min(shard_count);
    if shard_plans.is_empty() {
        return 1;
    }

    // Give each requested worker one partition-sized memory permit, while always
    // allowing the largest planned shard to make progress. Counting the largest
    // shards first makes the resulting concurrency safe for every scheduling order.
    let mut estimated_bytes: Vec<usize> = shard_plans
        .iter()
        .map(|plan| estimate_partitioned_shard_memory_bytes(plan.estimated_minimizers))
        .collect();
    estimated_bytes.sort_unstable_by(|left, right| right.cmp(left));
    let memory_budget: usize = requested_threads
        .saturating_mul(DEFAULT_PARTITION_TARGET_BYTES)
        .max(estimated_bytes[0]);
    let mut permitted: usize = 0;
    let mut reserved_bytes: usize = 0;
    for bytes in estimated_bytes.into_iter().take(thread_limit) {
        let Some(next_reserved) = reserved_bytes.checked_add(bytes) else {
            break;
        };
        if permitted > 0 && next_reserved > memory_budget {
            break;
        }
        reserved_bytes = next_reserved;
        permitted += 1;
    }

    permitted.max(1)
}

fn ceil_div_usize(numerator: usize, denominator: usize) -> usize {
    if denominator == 0 {
        return numerator;
    }

    numerator.saturating_add(denominator.saturating_sub(1)) / denominator
}

pub(crate) fn partition_build_plan(estimated_minimizers: usize) -> PartitionBuildPlan {
    let estimated_record_bytes: usize =
        estimated_minimizers.saturating_mul(size_of::<PartitionHitRecord>());
    let target_partition_bytes: usize =
        DEFAULT_PARTITION_TARGET_BYTES.max(size_of::<PartitionHitRecord>());
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

#[allow(dead_code)]
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
    let mut minimizer_window_count = 0;

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
    max_shard_minimizers: usize,
) -> io::Result<Vec<ShardPlan>> {
    validate_max_shard_minimizers(max_shard_minimizers)?;

    let mut plans = Vec::new();
    let mut shard_first_reference = 0;
    let mut shard_reference_count = 0;
    let mut shard_estimated_minimizers: usize = 0;
    for (reference_index, reference_minimizers) in minimizer_counts.iter().copied().enumerate() {
        let would_exceed_minimizers = shard_reference_count > 0
            && shard_estimated_minimizers.saturating_add(reference_minimizers)
                > max_shard_minimizers;

        if would_exceed_minimizers {
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
    max_shard_minimizers: usize,
    threads: usize,
    runtime_options: RuntimeOptions,
) -> io::Result<Vec<ShardPlan>> {
    validate_max_shard_minimizers(max_shard_minimizers)?;

    let plan_start: Instant = Instant::now();
    let planner_threads = threads.max(1).min(references.len().max(1));

    if runtime_options.progress_enabled {
        emit_progress(
            "shard_plan",
            &format!(
                "event=start\testimator=window-count\treferences={}\tmax_shard_minimizers={max_shard_minimizers}\tthreads={threads}\tplanner_parallelism={planner_threads}",
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
                Ok(reference_minimizers)
            })
            .collect::<io::Result<Vec<usize>>>()
    })?;

    let plans: Vec<ShardPlan> =
        plan_shards_from_minimizer_counts(&minimizer_counts, max_shard_minimizers)?;

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
    minimizer_hash_seed: u32,
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
        || manifest.minimizer_hash_seed != minimizer_hash_seed
        || manifest.key_mode != SKETCH_KEY_MODE
        || manifest.fragment_length != fragment_length
        || manifest.min_fragment_length != min_fragment_length
        || manifest.split_n_run != split_n_run
    {
        return Some(format!(
            "reference sketch database is incompatible: sketch_format_version={} database_schema_version={} k={} w={} minimizer_hash_seed={} key_mode={} fragment_length={} min_fragment_length={} split_n_run={}; rebuild the sketch with the requested parameters",
            manifest.sketch_format_version,
            manifest.database_schema_version,
            manifest.k,
            manifest.w,
            manifest.minimizer_hash_seed,
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
    partition_count: usize,
    /// Precomputed right-shift amount: `MinimizerKey::BITS - partition_count.trailing_zeros()`.
    /// Avoids recomputing `trailing_zeros()` for every minimizer pushed.
    partition_shift: u32,
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

        let partition_shift = if partition_count <= 1 {
            0u32
        } else {
            MinimizerKey::BITS - partition_count.trailing_zeros()
        };

        Ok(Self {
            files,
            writers,
            buffers: vec![Vec::new(); partition_count],
            buffer_record_limit: buffer_record_limit.max(1),
            partition_count,
            partition_shift,
        })
    }

    pub(crate) fn partition_count(&self) -> usize {
        self.files.len()
    }

    pub(crate) fn path(&self, partition_index: usize) -> &Path {
        &self.files[partition_index].path
    }

    #[allow(dead_code)]
    pub(crate) fn push(&mut self, record: PartitionHitRecord) -> io::Result<()> {
        let partition_index: usize = if self.partition_count <= 1 {
            0
        } else {
            (record.key >> self.partition_shift) as usize
        };
        let buffer: &mut Vec<PartitionHitRecord> = &mut self.buffers[partition_index];
        buffer.push(record);

        if buffer.len() >= self.buffer_record_limit {
            self.flush_partition(partition_index)?;
        }

        Ok(())
    }

    /// Push a slice of records using batch SIMD partition ID computation.
    /// Each partition is flushed as soon as it reaches `buffer_record_limit`, so a
    /// large input batch cannot cause an unbounded transient allocation.
    /// Re-benchmark this path in isolation before retuning or discarding it:
    /// a prior combined landing was masked by branch-predictor interference from
    /// unrelated hot-path changes rather than by a problem in batch routing itself.
    pub(crate) fn push_batch(&mut self, records: &[PartitionHitRecord]) -> io::Result<()> {
        if self.partition_count <= 1 {
            for chunk in records.chunks(self.buffer_record_limit) {
                if self.buffers[0].len().saturating_add(chunk.len()) > self.buffer_record_limit {
                    self.flush_partition(0)?;
                }
                self.buffers[0].extend_from_slice(chunk);
                if self.buffers[0].len() == self.buffer_record_limit {
                    self.flush_partition(0)?;
                }
            }
            return Ok(());
        }

        let shift = self.partition_shift;
        let limit = self.buffer_record_limit;
        let mut ids = [0usize; 8];
        let mut i = 0usize;

        // SIMD path: process 8 records at a time.
        while i + 8 <= records.len() {
            let chunk: &[PartitionHitRecord; 8] = records[i..i + 8]
                .try_into()
                .expect("slice of exactly 8 elements");
            compute_8_partition_ids(chunk, shift, &mut ids);
            for k in 0..8 {
                let partition_index = ids[k];
                self.buffers[partition_index].push(records[i + k]);
                if self.buffers[partition_index].len() == limit {
                    self.flush_partition(partition_index)?;
                }
            }
            i += 8;
        }

        // Scalar tail.
        for record in &records[i..] {
            let partition_index = (record.key >> shift) as usize;
            self.buffers[partition_index].push(*record);
            if self.buffers[partition_index].len() == limit {
                self.flush_partition(partition_index)?;
            }
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

/// Compute 8 partition IDs from 8 `PartitionHitRecord` keys simultaneously.
///
/// `shift` = `MinimizerKey::BITS - partition_bits`, precomputed in `PartitionWriters::new`.
///
/// `PartitionHitRecord` is `{ key: u32, hit: SeedHit { reference_contig_id: u32, position: u32 } }`
/// = 12 bytes, so keys are stride-3 u32s. They must be manually extracted before loading
/// into a SIMD register; a gather instruction handles this efficiently on AVX2 targets.
#[inline]
fn compute_8_partition_ids(records: &[PartitionHitRecord; 8], shift: u32, out: &mut [usize; 8]) {
    #[cfg(target_arch = "x86_64")]
    if std::arch::is_x86_feature_detected!("avx2") {
        // SAFETY: AVX2 confirmed above.
        unsafe { compute_8_ids_avx2(records, shift, out) };
        return;
    }

    // Scalar fallback, also the path on non-x86 targets.
    for (k, record) in records.iter().enumerate() {
        out[k] = (record.key >> shift) as usize;
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn compute_8_ids_avx2(records: &[PartitionHitRecord; 8], shift: u32, out: &mut [usize; 8]) {
    use core::arch::x86_64::{
        __m256i, _mm256_i32gather_epi32, _mm256_set1_epi32, _mm256_set_epi32, _mm256_srlv_epi32,
        _mm256_storeu_si256,
    };

    // PartitionHitRecord is 12 bytes = 3 x u32. Keys are at byte offsets
    // 0, 12, 24, 36, 48, 60, 72, 84 from the base of `records`.
    // _mm256_i32gather_epi32 with scale=1 uses raw byte offsets.
    let vindex = _mm256_set_epi32(84, 72, 60, 48, 36, 24, 12, 0);
    let base = records.as_ptr().cast::<i32>();

    // Gather: reads records[i].key for i in 0..8 in one instruction.
    let keys = unsafe { _mm256_i32gather_epi32::<1>(base, vindex) };

    // Logical right shift by `shift` bits to extract the partition index.
    // _mm256_srli_epi32 takes a compile-time-constant immediate in the intrinsic
    // form; use _mm256_srlv_epi32 for this runtime shift.
    let shift_vec = _mm256_set1_epi32(shift as i32);
    let shifted = _mm256_srlv_epi32(keys, shift_vec);

    // Store 8 x u32 results.
    let mut result = [0u32; 8];
    unsafe {
        _mm256_storeu_si256(result.as_mut_ptr().cast::<__m256i>(), shifted);
    }

    for k in 0..8 {
        out[k] = result[k] as usize;
    }
}

#[cfg(test)]
mod tests {
    use crate::ani::{
        database_build_parallelism, partition_build_plan, partition_id_for_key,
        plan_shards_from_minimizer_counts, PartitionBuildPlan, PartitionHitRecord,
        PartitionWriters, SeedHit, ShardPlan, MAX_PARTITION_COUNT, MIN_PARTITION_COUNT,
    };
    use std::io;

    #[test]
    fn push_batch_routes_identical_to_sequential_push() -> io::Result<()> {
        let partition_count = 16usize;
        let buffer_limit = 7;
        let tmp = std::env::temp_dir();
        let mut sequential = PartitionWriters::new(partition_count, Some(&tmp), buffer_limit)?;
        let mut batched = PartitionWriters::new(partition_count, Some(&tmp), buffer_limit)?;

        let records: Vec<PartitionHitRecord> = (0u32..200)
            .map(|i| PartitionHitRecord {
                key: i.wrapping_mul(0x9E37_79B9),
                hit: SeedHit {
                    reference_contig_id: i % 4,
                    position: i * 7,
                },
            })
            .collect();

        for record in &records {
            sequential.push(*record)?;
        }
        batched.push_batch(&records)?;

        for pid in 0..partition_count {
            assert_eq!(
                sequential.buffers[pid], batched.buffers[pid],
                "partition {pid} differs"
            );
            assert!(batched.buffers[pid].len() < buffer_limit);
        }

        Ok(())
    }

    #[test]
    fn partition_shift_is_precomputed_correctly() {
        for &partition_count in &[16usize, 64, 256, 1024, 4096] {
            let tmp = std::env::temp_dir();
            let pw = PartitionWriters::new(partition_count, Some(&tmp), 1).unwrap();
            let shift = pw.partition_shift;
            for key in [0u32, 1, u32::MAX / 2, u32::MAX - 1, u32::MAX] {
                let expected = partition_id_for_key(key, partition_count);
                let actual = (key >> shift) as usize;
                assert_eq!(actual, expected, "count={partition_count} key={key}");
            }
        }
    }

    #[test]
    fn shard_planner_splits_by_minimizer_target() -> io::Result<()> {
        let plans: Vec<ShardPlan> = plan_shards_from_minimizer_counts(&[200, 250, 100, 400], 500)?;

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
    fn shard_planner_allows_single_oversized_reference() -> io::Result<()> {
        let plans: Vec<ShardPlan> = plan_shards_from_minimizer_counts(&[600], 500)?;

        assert_eq!(
            plans,
            vec![ShardPlan {
                first_reference: 0,
                reference_count: 1,
                estimated_minimizers: 600,
            }]
        );

        Ok(())
    }

    #[test]
    fn automatic_partition_count_uses_power_of_two_bounds() {
        let small_plan: PartitionBuildPlan = partition_build_plan(10);
        assert_eq!(small_plan.partition_count, MIN_PARTITION_COUNT);
        assert!(small_plan.partition_count.is_power_of_two());

        let large_plan: PartitionBuildPlan = partition_build_plan(2_000_000_000);
        assert!(large_plan.partition_count > MIN_PARTITION_COUNT);
        assert!(large_plan.partition_count.is_power_of_two());
        assert!(large_plan.partition_count <= MAX_PARTITION_COUNT);
    }

    #[test]
    fn database_build_parallelism_uses_threads_when_memory_allows() {
        let shard_plans: Vec<ShardPlan> = vec![
            ShardPlan {
                first_reference: 0,
                reference_count: 1,
                estimated_minimizers: 1_000_000,
            };
            16
        ];
        assert_eq!(database_build_parallelism(12, &shard_plans), 12);
        assert_eq!(database_build_parallelism(4, &shard_plans), 4);
    }

    #[test]
    fn database_build_parallelism_limits_large_shards_by_memory_permits() {
        let shard_plans: Vec<ShardPlan> = vec![
            ShardPlan {
                first_reference: 0,
                reference_count: 1,
                estimated_minimizers: 300_000_000,
            };
            8
        ];

        assert_eq!(database_build_parallelism(8, &shard_plans), 1);
    }

    #[test]
    fn partition_id_uses_high_key_bits() {
        assert_eq!(partition_id_for_key(0x0000_0000, 16), 0);
        assert_eq!(partition_id_for_key(0x1000_0000, 16), 1);
        assert_eq!(partition_id_for_key(0xF000_0000, 16), 15);
    }
}
