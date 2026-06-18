//! Reference sketching, query mapping, and ANI summarization for `fasterANI`.
//!
//! The module keeps the FastANI-style data flow explicit:
//! reference minimizers are indexed once, query fragments are mapped through
//! seed-hit candidate regions and sliding-window identity scoring, and final output is summarized
//! per query/reference file pair.

use std::{
    cmp::Ordering,
    collections::{HashMap, HashSet},
    env, fs, io,
    io::{BufReader, BufWriter, Read, Write},
    mem::{align_of, size_of},
    os::fd::AsRawFd,
    path::{Path, PathBuf},
    ptr::NonNull,
    sync::{
        atomic::{AtomicUsize, Ordering as AtomicOrdering},
        Arc,
    },
    time::{Instant, SystemTime, UNIX_EPOCH},
};

use boomphf::Mphf;
use flate2::read::MultiGzDecoder;
use gzp::{deflate::Bgzf, ZBuilder};
use noodles::fasta;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use seq_hash::NtHasher;
use serde::{Deserialize, Serialize};
use simd_minimizers::canonical_minimizers;
use simd_minimizers::packed_seq::{PackedNSeqVec, Seq};

const DEFAULT_KMER_SIZE: usize = 16;
const DEFAULT_WINDOW_SIZE: usize = 24;
const MINIMIZER_HASH_SEED: u32 = 42;
const DEFAULT_FRAGMENT_LENGTH: u32 = 3000;
const DEFAULT_FRAGMENT_STRIDE: u32 = DEFAULT_FRAGMENT_LENGTH;
const DEFAULT_MIN_FRAGMENT_LENGTH: u32 = DEFAULT_FRAGMENT_LENGTH;
const DEFAULT_SPLIT_N_RUN: usize = 0;
const DEFAULT_MIN_PERCENT_IDENTITY: f64 = 80.0;
const DEFAULT_MASH_CONFIDENCE: f64 = 0.9;
const DEFAULT_FREQ_THRESHOLD_PERCENT: f64 = 0.0;
const DEFAULT_SHARD_SIZE: usize = 10_000;
const DEFAULT_SHARD_MINIMIZERS: usize = 500_000_000;
const PARTITIONED_INDEX_MINIMIZER_THRESHOLD: usize = 100_000_000;
const DEFAULT_PARTITION_TARGET_BYTES: usize = 512 * 1024 * 1024;
const MIN_PARTITION_COUNT: usize = 16;
const MAX_PARTITION_COUNT: usize = 4096;
const PARTITION_BUFFER_RECORDS: usize = 65_536;
const ESTIMATED_PARTITIONED_SHARD_BYTES_PER_MINIMIZER: usize = 16;
const SKETCH_MAGIC: &[u8; 8] = b"FANIIDX1";
const SKETCH_VERSION: u32 = 9;
const SKETCH_DATABASE_SCHEMA_VERSION: u32 = 1;
const SKETCH_KEY_MODE: &str = "canonical-2bit-u32";
const REFERENCE_PROGRESS_INTERVAL: usize = 1000;
const SKETCH_KEY_PACK_PROGRESS_INTERVAL: usize = 5_000_000;
#[cfg(test)]
const SKETCH_CONTIG_PACK_PROGRESS_INTERVAL: usize = 100_000;

type MinimizerKey = u32;
type ReferenceHitMap = FxHashMap<MinimizerKey, Vec<SeedHit>>;

/// User-selectable strategy for building the reference sketch index.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum IndexBuildMode {
    Auto,
    Hash,
    Partitioned,
}

impl IndexBuildMode {
    fn name(self) -> &'static str {
        match self {
            Self::Auto => "auto",
            Self::Hash => "hash",
            Self::Partitioned => "partitioned",
        }
    }

    fn parse(value: &str) -> io::Result<Self> {
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
struct PartitionHitRecord {
    key: MinimizerKey,
    hit: SeedHit,
}

/// One grouped minimizer key and the hit range assigned to it in the final payload file.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct GroupedKeyRecord {
    key: MinimizerKey,
    hit_offset: u64,
    hit_count: u32,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct PartitionBuildPlan {
    partition_count: usize,
    target_partition_bytes: usize,
    estimated_record_bytes: usize,
}

struct PartitionGroupResult {
    partition_index: usize,
    grouped_key_scratch: ScratchFile,
    hit_payload_scratch: ScratchFile,
    keys: Vec<MinimizerKey>,
    hit_count: usize,
}

/// One input reference genome file and the number of bases FastANI considers mappable.
#[derive(Clone, Deserialize, Serialize)]
struct ReferenceFile {
    path: String,
    mapped_length: u64,
}

/// One FASTA record from a reference file and its ordered minimizers.
#[derive(Clone, Deserialize, Serialize)]
struct ReferenceContig {
    file_id: usize,
    minimizers: Vec<ReferenceMinimizer>,
}

/// Optional human-readable metadata for one indexed reference contig or split segment.
#[derive(Clone)]
struct ReferenceContigName {
    file_id: usize,
    name: String,
    segment_start: u32,
    segment_end: u32,
}

/// Compact persisted metadata for one reference contig in a memory-mapped sketch.
#[repr(C)]
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct ContigRecord {
    minimizer_offset: u64,
    file_id: u32,
    minimizer_count: u32,
}

/// A reference minimizer hash and its zero-based contig position.
#[repr(C)]
#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
struct ReferenceMinimizer {
    hash: MinimizerKey,
    position: u32,
}

/// Compact seed hit stored in the reference index.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, Deserialize, Eq, PartialEq, Serialize)]
struct SeedHit {
    reference_contig_id: u32,
    position: u32,
}

/// Fully loaded reference sketch, backed by either an in-memory hash map or an mmap cache.
struct ReferenceSketch {
    files: Vec<ReferenceFile>,
    contigs: ReferenceContigs,
    contig_names: Option<Vec<ReferenceContigName>>,
    index: ReferenceIndex,
}

/// Reference contig minimizers, either owned after a fresh build or sliced from a loaded sketch.
enum ReferenceContigs {
    Owned(Vec<ReferenceContig>),
    Mmap(MmapReferenceContigs),
}

/// Lookup table from minimizer hash to all reference positions containing that minimizer.
enum ReferenceIndex {
    Hash(ReferenceHitMap),
    Mphf(MmapReferenceIndex),
}

/// JSON metadata stored at the front of a `.fasketch` cache.
#[derive(Deserialize, Serialize)]
struct CachedReferenceMetadata {
    version: u32,
    k: usize,
    w: usize,
    #[serde(default)]
    key_mode: String,
    #[serde(default = "default_fragment_length")]
    fragment_length: u32,
    min_fragment_length: u32,
    #[serde(default)]
    split_n_run: usize,
    #[serde(default)]
    dust_enabled: bool,
    files: Vec<ReferenceFile>,
    mphf: Mphf<MinimizerKey>,
    key_count: usize,
    hit_count: usize,
    contig_count: usize,
    reference_minimizer_count: usize,
}

/// Top-level metadata for a manifest-backed sharded reference database.
#[derive(Clone, Deserialize, Serialize)]
struct ShardManifest {
    sketch_format_version: u32,
    database_schema_version: u32,
    k: usize,
    w: usize,
    key_mode: String,
    #[serde(default = "default_fragment_length")]
    fragment_length: u32,
    min_fragment_length: u32,
    split_n_run: usize,
    #[serde(default)]
    dust_enabled: bool,
    shard_size: usize,
    #[serde(default = "default_shard_minimizers")]
    shard_minimizers: usize,
    total_references: usize,
    total_reference_contigs: usize,
    total_mapped_reference_length: u64,
    total_reference_minimizers: usize,
    total_shard_unique_minimizers: usize,
    build_unix_seconds: u64,
    build_args: Vec<String>,
    reference_list_checksum: u64,
    shards: Vec<ShardManifestEntry>,
}

/// One physical sketch shard referenced by a sharded database manifest.
#[derive(Clone, Deserialize, Serialize)]
struct ShardManifestEntry {
    shard_index: usize,
    filename: String,
    first_reference: usize,
    reference_count: usize,
    reference_contigs: usize,
    mapped_reference_length: u64,
    reference_minimizers: usize,
    unique_minimizers: usize,
}

/// Small summary returned by a streaming shard build.
#[derive(Clone, Copy, Default)]
struct SketchBuildStats {
    reference_count: usize,
    reference_contig_count: usize,
    mapped_reference_length: u64,
    reference_minimizer_count: usize,
    unique_minimizer_count: usize,
}

/// Completed shard build result used to assemble the final database manifest.
struct ShardBuildResult {
    entry: ShardManifestEntry,
}

/// Reference-list range selected for one physical sketch shard.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct ShardPlan {
    first_reference: usize,
    reference_count: usize,
    estimated_minimizers: usize,
}

/// Reference database opened by the CLI, either legacy single-sketch or manifest-backed shards.
enum SketchDatabase {
    Single(ReferenceSketch),
    Sharded {
        prefix: PathBuf,
        manifest: ShardManifest,
    },
}

#[cfg(debug_assertions)]
#[derive(Default)]
struct ReferenceMemoryEstimate {
    reference_minimizers: usize,
    reference_minimizer_vec_bytes: usize,
    unique_index_keys: usize,
    seed_hits: usize,
    seed_hit_vec_bytes: usize,
    hash_index_rough_bytes: usize,
    mmap_file_bytes: usize,
    mmap_slot_key_bytes: usize,
    mmap_hit_offset_bytes: usize,
    mmap_hit_count_bytes: usize,
    mmap_hit_payload_bytes: usize,
    mmap_contig_record_bytes: usize,
    mmap_reference_minimizer_bytes: usize,
}

/// Runtime controls shared by long-running reference build and sketch operations.
#[derive(Clone, Copy, Default)]
struct RuntimeOptions {
    progress_enabled: bool,
    max_memory_bytes: Option<u64>,
    worker_threads: usize,
}

impl RuntimeOptions {
    fn effective_worker_threads(self) -> usize {
        self.worker_threads.max(1)
    }

    fn with_worker_threads(self, worker_threads: usize) -> Self {
        Self {
            worker_threads: worker_threads.max(1),
            ..self
        }
    }
}

#[cfg(debug_assertions)]
#[derive(Default)]
struct QueryMemoryEstimate {
    fragment_struct_bytes: usize,
    query_minimizer_vec_bytes: usize,
    seed_minimizer_vec_bytes: usize,
}

impl ReferenceIndex {
    fn get(&self, minimizer: &MinimizerKey) -> Option<&[SeedHit]> {
        match self {
            Self::Hash(index) => index.get(minimizer).map(Vec::as_slice),
            Self::Mphf(index) => index.get(minimizer),
        }
    }

    fn len(&self) -> usize {
        match self {
            Self::Hash(index) => index.len(),
            Self::Mphf(index) => index.key_count,
        }
    }

    fn frequency_threshold(&self, percent: f64) -> usize {
        if percent <= 0.0 {
            return usize::MAX;
        }

        let mut histogram: HashMap<usize, usize> = HashMap::new();
        let mut total_unique_minimizers: usize = 0usize;

        match self {
            Self::Hash(index) => {
                for hits in index.values() {
                    *histogram.entry(hits.len()).or_default() += 1;
                    total_unique_minimizers += 1;
                }
            }
            Self::Mphf(index) => {
                for &count in index.hit_counts() {
                    *histogram.entry(count as usize).or_default() += 1;
                    total_unique_minimizers += 1;
                }
            }
        }

        let minimizers_to_ignore: usize =
            (total_unique_minimizers as f64 * percent / 100.0) as usize;
        if minimizers_to_ignore == 0 {
            return usize::MAX;
        }

        let mut frequencies: Vec<(usize, usize)> = histogram.into_iter().collect::<Vec<_>>();
        frequencies.sort_unstable_by(|left, right| right.0.cmp(&left.0));

        let mut sum: usize = 0usize;
        let mut threshold: usize = usize::MAX;
        for (frequency, count) in frequencies {
            sum += count;
            if sum < minimizers_to_ignore {
                threshold = frequency;
            } else if sum == minimizers_to_ignore {
                threshold = frequency;
                break;
            } else {
                break;
            }
        }

        threshold
    }
}

/// Memory-mapped minimal-perfect-hash index and the arrays it addresses.
struct MmapReferenceIndex {
    mphf: Mphf<MinimizerKey>,
    mmap: Arc<MmapFile>,
    key_count: usize,
    hit_count: usize,
    slot_keys_offset: usize,
    hit_offsets_offset: usize,
    hit_counts_offset: usize,
    hit_payloads_offset: usize,
}

impl MmapReferenceIndex {
    fn get(&self, minimizer: &MinimizerKey) -> Option<&[SeedHit]> {
        let slot = self.mphf.try_hash(minimizer)? as usize;
        if slot >= self.key_count {
            return None;
        }

        let slot_keys = self.slot_keys();
        if slot_keys[slot] != *minimizer {
            return None;
        }

        let start = self.hit_offsets()[slot] as usize;
        let count = self.hit_counts()[slot] as usize;
        let end = start.checked_add(count)?;
        self.hit_payloads().get(start..end)
    }

    fn slot_keys(&self) -> &[MinimizerKey] {
        mmap_slice_at(&self.mmap, self.slot_keys_offset, self.key_count)
    }

    fn hit_offsets(&self) -> &[u64] {
        mmap_slice_at(&self.mmap, self.hit_offsets_offset, self.key_count)
    }

    fn hit_counts(&self) -> &[u32] {
        mmap_slice_at(&self.mmap, self.hit_counts_offset, self.key_count)
    }

    fn hit_payloads(&self) -> &[SeedHit] {
        mmap_slice_at(&self.mmap, self.hit_payloads_offset, self.hit_count)
    }
}

/// Memory-mapped contig records and the flat minimizer array they address.
struct MmapReferenceContigs {
    mmap: Arc<MmapFile>,
    contig_count: usize,
    reference_minimizer_count: usize,
    contig_records_offset: usize,
    reference_minimizers_offset: usize,
}

impl MmapReferenceContigs {
    fn records(&self) -> &[ContigRecord] {
        mmap_slice_at(&self.mmap, self.contig_records_offset, self.contig_count)
    }

    fn reference_minimizers(&self) -> &[ReferenceMinimizer] {
        mmap_slice_at(
            &self.mmap,
            self.reference_minimizers_offset,
            self.reference_minimizer_count,
        )
    }
}

impl ReferenceContigs {
    fn len(&self) -> usize {
        match self {
            Self::Owned(contigs) => contigs.len(),
            Self::Mmap(contigs) => contigs.contig_count,
        }
    }

    fn file_id(&self, contig_id: usize) -> Option<usize> {
        match self {
            Self::Owned(contigs) => contigs.get(contig_id).map(|contig| contig.file_id),
            Self::Mmap(contigs) => contigs
                .records()
                .get(contig_id)
                .map(|record| record.file_id as usize),
        }
    }

    fn minimizers(&self, contig_id: usize) -> Option<&[ReferenceMinimizer]> {
        match self {
            Self::Owned(contigs) => contigs
                .get(contig_id)
                .map(|contig| contig.minimizers.as_slice()),
            Self::Mmap(contigs) => {
                let record: ContigRecord = *contigs.records().get(contig_id)?;
                let start: usize = record.minimizer_offset as usize;
                let count: usize = record.minimizer_count as usize;
                let end: usize = start.checked_add(count)?;
                contigs.reference_minimizers().get(start..end)
            }
        }
    }

    #[cfg(debug_assertions)]
    fn total_minimizers(&self) -> usize {
        match self {
            Self::Owned(contigs) => contigs.iter().map(|contig| contig.minimizers.len()).sum(),
            Self::Mmap(contigs) => contigs.reference_minimizer_count,
        }
    }

    #[cfg(debug_assertions)]
    fn owned_minimizer_capacity_bytes(&self) -> usize {
        match self {
            Self::Owned(contigs) => contigs
                .iter()
                .map(|contig| contig.minimizers.capacity() * size_of::<ReferenceMinimizer>())
                .sum(),
            Self::Mmap(_) => 0,
        }
    }
}

fn mmap_slice_at<T>(mmap: &MmapFile, offset: usize, count: usize) -> &[T] {
    let byte_len: usize = count
        .checked_mul(size_of::<T>())
        .expect("mmap slice length overflow");
    let end: usize = offset
        .checked_add(byte_len)
        .expect("mmap slice offset overflow");
    assert!(end <= mmap.as_slice().len());
    assert_eq!(
        (mmap.as_slice().as_ptr() as usize + offset) % align_of::<T>(),
        0
    );

    unsafe { std::slice::from_raw_parts(mmap.as_slice().as_ptr().add(offset) as *const T, count) }
}

/// Read-only memory map wrapper for cached reference sketches.
struct MmapFile {
    ptr: NonNull<u8>,
    len: usize,
}

unsafe impl Send for MmapFile {}
unsafe impl Sync for MmapFile {}

impl MmapFile {
    fn open(path: &Path) -> io::Result<Self> {
        let file: fs::File = fs::File::open(path)?;
        let len: usize = file.metadata()?.len() as usize;
        if len == 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "cannot mmap an empty reference sketch",
            ));
        }

        let ptr: *mut libc::c_void = unsafe {
            libc::mmap(
                std::ptr::null_mut(),
                len,
                libc::PROT_READ,
                libc::MAP_PRIVATE,
                file.as_raw_fd(),
                0,
            )
        };

        if ptr == libc::MAP_FAILED {
            return Err(io::Error::last_os_error());
        }

        Ok(Self {
            ptr: NonNull::new(ptr as *mut u8).expect("mmap returned null"),
            len,
        })
    }

    fn as_slice(&self) -> &[u8] {
        unsafe { std::slice::from_raw_parts(self.ptr.as_ptr(), self.len) }
    }
}

impl Drop for MmapFile {
    fn drop(&mut self) {
        unsafe {
            libc::munmap(self.ptr.as_ptr() as *mut libc::c_void, self.len);
        }
    }
}

/// One fixed-length query fragment with full scoring minimizers and optional reduced seed hashes.
#[derive(Clone)]
struct QueryFragment {
    id: usize,
    contig_id: usize,
    start: u32,
    end: u32,
    length: u32,
    minimizers: Vec<MinimizerKey>,
    seed_minimizers: Vec<MinimizerKey>,
}

/// Query genome split into FastANI-style fixed-width fragments.
struct QueryFile {
    fragments: Vec<QueryFragment>,
    contig_names: Vec<String>,
    mapped_length: u64,
}

/// Candidate reference position interval produced by clustered seed hits.
#[derive(Clone, Copy)]
struct ReferenceCandidateRegion {
    reference_contig_id: usize,
    start_position: u32,
    end_position: u32,
}

/// One retained mapping between a query fragment and a reference location.
#[derive(Clone)]
struct MappingResult {
    reference_file_id: usize,
    reference_contig_id: usize,
    query_fragment_id: usize,
    query_fragment_length: u32,
    reference_start: u32,
    identity: f64,
    query_minimizer_count: usize,
    reference_minimizer_count: usize,
    shared_minimizers: usize,
    union_minimizers: usize,
    jaccard: f64,
}

/// Stable identity for marking mappings that survive reciprocal-best filtering.
#[derive(Hash, Eq, PartialEq)]
struct MappingResultKey {
    reference_file_id: usize,
    reference_contig_id: usize,
    query_fragment_id: usize,
    query_fragment_length: u32,
    reference_start: u32,
    identity_bits: u64,
    query_minimizer_count: usize,
    reference_minimizer_count: usize,
    shared_minimizers: usize,
    union_minimizers: usize,
}

impl MappingResultKey {
    fn from_mapping(mapping: &MappingResult) -> Self {
        Self {
            reference_file_id: mapping.reference_file_id,
            reference_contig_id: mapping.reference_contig_id,
            query_fragment_id: mapping.query_fragment_id,
            query_fragment_length: mapping.query_fragment_length,
            reference_start: mapping.reference_start,
            identity_bits: mapping.identity.to_bits(),
            query_minimizer_count: mapping.query_minimizer_count,
            reference_minimizer_count: mapping.reference_minimizer_count,
            shared_minimizers: mapping.shared_minimizers,
            union_minimizers: mapping.union_minimizers,
        }
    }
}

/// Accumulator for final ANI output for one reference file.
#[derive(Default)]
struct AniSummary {
    shared_fragments: usize,
    shared_bases: u64,
    weighted_identity_sum: f64,
}

/// Final per-reference summaries plus the exact mappings that contributed to them.
struct AniComputation {
    summaries: Vec<AniSummary>,
    reciprocal_best_keys: HashSet<MappingResultKey>,
}

/// Per-thread reusable buffers for query-to-reference mapping.
#[derive(Default)]
struct MappingScratch {
    seed_hits: Vec<SeedHit>,
    candidate_regions: Vec<ReferenceCandidateRegion>,
    fragment_mappings: Vec<MappingResult>,
    counter: SlidingSketchCounter,
}

/// Fine-grained counters for the query mapping hot path.
#[cfg(debug_assertions)]
const SEED_HIT_HISTOGRAM_UPPER_BOUNDS: [usize; 13] =
    [1, 2, 4, 9, 24, 49, 99, 249, 499, 999, 4_999, 9_999, 49_999];

#[cfg(debug_assertions)]
const SEED_HIT_HISTOGRAM_OVERFLOW_LABEL: &str = "50000+";

#[cfg(debug_assertions)]
#[derive(Clone, Default)]
struct MappingMetrics {
    candidate_discovery_calls: usize,
    seed_hits_collected: usize,
    candidate_regions_found: usize,
    candidate_regions_scored: usize,
    reference_minimizers_scanned: usize,
    scoring_window_steps: usize,
    retained_mappings: usize,
    candidate_discovery_elapsed: std::time::Duration,
    scoring_elapsed: std::time::Duration,
    seed_lookup_count: usize,
    seed_lookup_zero_hits: usize,
    seed_lookup_skipped_by_frequency: usize,
    seed_hit_list_max: usize,
    seed_hit_list_bins: [usize; SEED_HIT_HISTOGRAM_UPPER_BOUNDS.len() + 1],
    seed_hit_list_bin_hits: [usize; SEED_HIT_HISTOGRAM_UPPER_BOUNDS.len() + 1],
}

#[cfg(debug_assertions)]
impl MappingMetrics {
    fn merge(&mut self, other: Self) {
        self.candidate_discovery_calls += other.candidate_discovery_calls;
        self.seed_hits_collected += other.seed_hits_collected;
        self.candidate_regions_found += other.candidate_regions_found;
        self.candidate_regions_scored += other.candidate_regions_scored;
        self.reference_minimizers_scanned += other.reference_minimizers_scanned;
        self.scoring_window_steps += other.scoring_window_steps;
        self.retained_mappings += other.retained_mappings;
        self.candidate_discovery_elapsed += other.candidate_discovery_elapsed;
        self.scoring_elapsed += other.scoring_elapsed;
        self.seed_lookup_count += other.seed_lookup_count;
        self.seed_lookup_zero_hits += other.seed_lookup_zero_hits;
        self.seed_lookup_skipped_by_frequency += other.seed_lookup_skipped_by_frequency;
        self.seed_hit_list_max = self.seed_hit_list_max.max(other.seed_hit_list_max);

        for (left, right) in self
            .seed_hit_list_bins
            .iter_mut()
            .zip(other.seed_hit_list_bins)
        {
            *left += right;
        }

        for (left, right) in self
            .seed_hit_list_bin_hits
            .iter_mut()
            .zip(other.seed_hit_list_bin_hits)
        {
            *left += right;
        }
    }

    fn record_seed_lookup(&mut self, hit_list_len: Option<usize>, frequency_threshold: usize) {
        self.seed_lookup_count += 1;

        let Some(hit_list_len) = hit_list_len else {
            self.seed_lookup_zero_hits += 1;
            return;
        };

        if hit_list_len >= frequency_threshold {
            self.seed_lookup_skipped_by_frequency += 1;
        }

        self.seed_hit_list_max = self.seed_hit_list_max.max(hit_list_len);
        let bin_index: usize = seed_hit_histogram_bin_index(hit_list_len);
        self.seed_hit_list_bins[bin_index] += 1;
        self.seed_hit_list_bin_hits[bin_index] += hit_list_len;
    }
}

#[cfg(not(debug_assertions))]
#[derive(Clone, Default)]
struct MappingMetrics;

#[cfg(debug_assertions)]
fn seed_hit_histogram_bin_index(hit_list_len: usize) -> usize {
    SEED_HIT_HISTOGRAM_UPPER_BOUNDS
        .iter()
        .position(|&upper_bound| hit_list_len <= upper_bound)
        .unwrap_or(SEED_HIT_HISTOGRAM_UPPER_BOUNDS.len())
}

/// Raw mapping results plus optional hot-path metrics for one query file.
struct MappingOutput {
    results: Vec<MappingResult>,
    #[cfg(debug_assertions)]
    metrics: MappingMetrics,
}

/// Sliding minimizer set used while scoring candidate reference windows.
#[derive(Default)]
struct SlidingSketchCounter {
    hashes: Vec<MinimizerKey>,
    query_present: Vec<bool>,
    reference_counts: Vec<usize>,
    active: Fenwick,
    shared: Fenwick,
    sketch_size: usize,
}

/// Minimizers observed for one sequence plus the raw minimizer-window coverage signal.
struct MinimizerObservation {
    minimizers_with_positions: Vec<(MinimizerKey, u32)>,
}

/// Query-fragment sketch with both scoring minimizers and the missing-window quality signal.
struct QueryFragmentSketch {
    minimizers: Vec<MinimizerKey>,
    seed_minimizers: Vec<MinimizerKey>,
}

/// Fenwick tree for rank/select operations over active minimizer coordinates.
#[derive(Default)]
struct Fenwick {
    tree: Vec<i32>,
    total: i32,
}

impl Fenwick {
    fn reset(&mut self, len: usize) {
        self.tree.clear();
        self.tree.resize(len + 1, 0);
        self.total = 0;
    }

    fn add(&mut self, index: usize, delta: i32) {
        let mut i = index + 1;
        self.total += delta;

        while i < self.tree.len() {
            self.tree[i] += delta;
            i += i & i.wrapping_neg();
        }
    }

    fn prefix_sum_inclusive(&self, index: usize) -> i32 {
        let mut i = index + 1;
        let mut sum = 0;

        while i > 0 {
            sum += self.tree[i];
            i &= i - 1;
        }

        sum
    }

    fn select_by_rank(&self, rank: i32) -> Option<usize> {
        if rank <= 0 || rank > self.total {
            return None;
        }

        let mut idx = 0usize;
        let mut bit = 1usize;

        while bit < self.tree.len() {
            bit <<= 1;
        }

        let mut remaining = rank;

        while bit > 0 {
            let next = idx + bit;

            if next < self.tree.len() && self.tree[next] < remaining {
                idx = next;
                remaining -= self.tree[next];
            }

            bit >>= 1;
        }

        Some(idx)
    }
}

impl SlidingSketchCounter {
    fn clear(&mut self) {
        for count in &mut self.reference_counts {
            *count = 0;
        }

        self.active.reset(self.hashes.len());
        self.shared.reset(self.hashes.len());

        for (index, &present) in self.query_present.iter().enumerate() {
            if present {
                self.active.add(index, 1);
            }
        }
    }

    fn prepare(
        &mut self,
        query_minimizers: &[MinimizerKey],
        reference_minimizers: &[ReferenceMinimizer],
    ) {
        self.hashes.clear();
        self.hashes.extend_from_slice(query_minimizers);
        self.hashes
            .extend(reference_minimizers.iter().map(|minimizer| minimizer.hash));
        self.hashes.sort_unstable();
        self.hashes.dedup();

        self.query_present.clear();
        self.query_present.resize(self.hashes.len(), false);

        for minimizer in query_minimizers {
            let index = self
                .hashes
                .binary_search(minimizer)
                .expect("query minimizer missing from coordinate table");
            self.query_present[index] = true;
        }

        self.reference_counts.clear();
        self.reference_counts.resize(self.hashes.len(), 0);
        self.sketch_size = query_minimizers.len();
        self.clear();
    }

    fn insert(&mut self, hash: MinimizerKey) {
        let index = self
            .hashes
            .binary_search(&hash)
            .expect("reference minimizer missing from coordinate table");
        let count = &mut self.reference_counts[index];

        if *count == 0 {
            if self.query_present[index] {
                self.shared.add(index, 1);
            } else {
                self.active.add(index, 1);
            }
        }

        *count += 1;
    }

    fn remove(&mut self, hash: MinimizerKey) {
        let index = self
            .hashes
            .binary_search(&hash)
            .expect("reference minimizer missing from coordinate table");
        let count = &mut self.reference_counts[index];

        if *count == 0 {
            return;
        }

        *count -= 1;

        if *count == 0 {
            if self.query_present[index] {
                self.shared.add(index, -1);
            } else {
                self.active.add(index, -1);
            }
        }
    }

    fn shared_count(&self) -> usize {
        let rank = self.sketch_size as i32;
        let Some(pivot) = self.active.select_by_rank(rank) else {
            return 0;
        };

        self.shared.prefix_sum_inclusive(pivot) as usize
    }

    fn reference_minimizer_count(&self) -> usize {
        self.reference_counts
            .iter()
            .filter(|&&count| count > 0)
            .count()
    }
}

fn expected_minimizer_window_count(
    sequence_len: usize,
    kmer_size: usize,
    window_size: usize,
) -> usize {
    let Some(minimum_sequence_len) = kmer_size.checked_add(window_size.saturating_sub(1)) else {
        return 0;
    };

    if sequence_len < minimum_sequence_len {
        0
    } else {
        sequence_len - minimum_sequence_len + 1
    }
}

#[cfg(test)]
fn usable_minimizer_window_count(sequence: &[u8], kmer_size: usize, window_size: usize) -> usize {
    let expected_window_count: usize =
        expected_minimizer_window_count(sequence.len(), kmer_size, window_size);
    if expected_window_count == 0 {
        return 0;
    }

    let minimizer_window_span: usize = kmer_size + window_size - 1;
    let mut ambiguous_prefix: Vec<usize> = Vec::with_capacity(sequence.len() + 1);
    ambiguous_prefix.push(0);

    for base in sequence {
        let previous_count: usize = *ambiguous_prefix
            .last()
            .expect("ambiguous prefix always has a zero entry");
        let next_count: usize = previous_count + usize::from(matches!(base, b'N' | b'n'));
        ambiguous_prefix.push(next_count);
    }

    (0..expected_window_count)
        .filter(|&start| ambiguous_prefix[start + minimizer_window_span] == ambiguous_prefix[start])
        .count()
}

/// Return canonical minimizers plus how many minimizer windows were usable before deduplication.
fn canonical_minimizer_observation(
    sequence: &[u8],
    kmer_size: usize,
    window_size: usize,
) -> MinimizerObservation {
    if sequence.len() < kmer_size + window_size.saturating_sub(1) || is_all_n_sequence(sequence) {
        return MinimizerObservation {
            minimizers_with_positions: Vec::new(),
        };
    }

    let packed_sequence: PackedNSeqVec = PackedNSeqVec::from_ascii(sequence);
    let packed_sequence_slice = packed_sequence.as_slice();
    let sequence_slice = packed_sequence_slice.seq;
    let hasher: NtHasher<true> = NtHasher::<true>::new_with_seed(kmer_size, MINIMIZER_HASH_SEED);
    let mut minimizer_positions: Vec<u32> = Vec::new();
    let minimizer_builder = canonical_minimizers(kmer_size, window_size).hasher(&hasher);
    let _ = minimizer_builder
        .run_skip_ambiguous_windows(packed_sequence_slice, &mut minimizer_positions);
    debug_assert!(kmer_size <= 16, "u32 2-bit minimizer keys require k <= 16");

    let minimizers_with_positions: Vec<(MinimizerKey, u32)> = minimizer_positions
        .into_iter()
        .filter_map(|position| {
            let pos = position as usize;
            let forward_kmer = sequence_slice.read_kmer(kmer_size, pos);
            let reverse_complement_kmer = sequence_slice.read_revcomp_kmer(kmer_size, pos);
            if forward_kmer == reverse_complement_kmer {
                return None;
            }

            let canonical_kmer: u64 = forward_kmer.min(reverse_complement_kmer);
            let key: MinimizerKey = MinimizerKey::try_from(canonical_kmer)
                .expect("u32 2-bit minimizer keys require k <= 16");

            Some((key, position))
        })
        .collect();

    MinimizerObservation {
        minimizers_with_positions,
    }
}

/// Return canonical minimizer hashes and their positions for one nucleotide sequence.
fn canonical_minimizers_with_positions(
    sequence: &[u8],
    kmer_size: usize,
    window_size: usize,
) -> Vec<(MinimizerKey, u32)> {
    canonical_minimizer_observation(sequence, kmer_size, window_size).minimizers_with_positions
}

fn is_all_n_sequence(sequence: &[u8]) -> bool {
    !sequence.is_empty() && sequence.iter().all(|base| matches!(base, b'N' | b'n'))
}

fn is_ambiguous_base(base: u8) -> bool {
    matches!(base, b'N' | b'n')
}

fn is_no_usable_fragments_error(error: &io::Error) -> bool {
    error.kind() == io::ErrorKind::InvalidData
        && error.to_string() == "ERROR: Input has no usable fragments"
}

/// Return the minimizer hashes used for candidate discovery.
fn select_seed_minimizers(
    minimizers: &[MinimizerKey],
    minmer_count: Option<usize>,
) -> Vec<MinimizerKey> {
    match minmer_count {
        Some(count) => minimizers.iter().take(count).copied().collect(),
        None => minimizers.to_vec(),
    }
}

fn query_fragment_sketch(
    fragment_sequence: &[u8],
    kmer_size: usize,
    window_size: usize,
    minmer_count: Option<usize>,
) -> QueryFragmentSketch {
    let observation: MinimizerObservation =
        canonical_minimizer_observation(fragment_sequence, kmer_size, window_size);
    let mut minimizers: Vec<MinimizerKey> = observation
        .minimizers_with_positions
        .into_iter()
        .map(|(minimizer, _position)| minimizer)
        .collect();
    minimizers.sort_unstable();
    minimizers.dedup();
    let seed_minimizers: Vec<MinimizerKey> = select_seed_minimizers(&minimizers, minmer_count);

    QueryFragmentSketch {
        minimizers,
        seed_minimizers,
    }
}

fn split_sequence_ranges(sequence: &[u8], split_n_run: usize) -> Vec<std::ops::Range<usize>> {
    if sequence.is_empty() {
        return Vec::new();
    }

    if split_n_run == 0 {
        return vec![0..sequence.len()];
    }

    let mut ranges: Vec<std::ops::Range<usize>> = Vec::new();
    let mut segment_start: usize = 0usize;
    let mut position: usize = 0usize;

    while position < sequence.len() {
        if !is_ambiguous_base(sequence[position]) {
            position += 1;
            continue;
        }

        let run_start: usize = position;
        while position < sequence.len() && is_ambiguous_base(sequence[position]) {
            position += 1;
        }
        let run_end: usize = position;

        if run_end - run_start >= split_n_run {
            if segment_start < run_start {
                ranges.push(segment_start..run_start);
            }
            segment_start = run_end;
        }
    }

    if segment_start < sequence.len() {
        ranges.push(segment_start..sequence.len());
    }

    ranges
}

/// Return query fragment ranges for either FastANI-compatible chunks or adaptive overlap mode.
fn query_fragment_ranges(
    sequence_len: usize,
    fragment_length: usize,
    fragment_stride: usize,
    min_fragment_length: usize,
) -> Vec<std::ops::Range<usize>> {
    if sequence_len < min_fragment_length {
        return Vec::new();
    }

    let retain_tail_fragment: bool = min_fragment_length < fragment_length;
    let fastani_compatible: bool = fragment_stride == fragment_length && !retain_tail_fragment;

    if fastani_compatible {
        return (0..sequence_len)
            .step_by(fragment_length)
            .take_while(|start| start + fragment_length <= sequence_len)
            .map(|start| start..start + fragment_length)
            .collect();
    }

    if sequence_len <= fragment_length {
        return vec![0..sequence_len];
    }

    let mut ranges: Vec<std::ops::Range<usize>> = Vec::new();
    let mut start: usize = 0usize;

    while start + fragment_length <= sequence_len {
        ranges.push(start..start + fragment_length);
        start += fragment_stride;
    }

    let tail_len: usize = sequence_len.saturating_sub(start);
    if tail_len >= min_fragment_length {
        ranges.push(start..sequence_len);
        return ranges;
    }

    if retain_tail_fragment {
        return ranges;
    }

    let tail_start: usize = sequence_len - fragment_length;
    if ranges
        .last()
        .is_none_or(|last_range| last_range.start != tail_start)
    {
        ranges.push(tail_start..sequence_len);
    }

    ranges
}

fn mapped_length_from_fragment_ranges(
    sequence_len: usize,
    fragment_length: u32,
    min_fragment_length: u32,
) -> u64 {
    query_fragment_ranges(
        sequence_len,
        fragment_length as usize,
        fragment_length as usize,
        min_fragment_length as usize,
    )
    .into_iter()
    .map(|range| range.len() as u64)
    .sum()
}

fn fastani_compatible_fragment_mode(
    fragment_length: u32,
    fragment_stride: u32,
    min_fragment_length: u32,
) -> bool {
    fragment_stride == fragment_length && min_fragment_length == fragment_length
}

#[cfg(test)]
mod tests {
    use super::*;

    fn repeated_acgt(len: usize) -> Vec<u8> {
        (0..len).map(|i| b"ACGT"[i % 4]).collect()
    }

    fn sample_shard_manifest() -> ShardManifest {
        ShardManifest {
            sketch_format_version: SKETCH_VERSION,
            database_schema_version: SKETCH_DATABASE_SCHEMA_VERSION,
            k: DEFAULT_KMER_SIZE,
            w: DEFAULT_WINDOW_SIZE,
            key_mode: SKETCH_KEY_MODE.to_string(),
            fragment_length: DEFAULT_FRAGMENT_LENGTH,
            min_fragment_length: DEFAULT_MIN_FRAGMENT_LENGTH,
            split_n_run: DEFAULT_SPLIT_N_RUN,
            dust_enabled: false,
            shard_size: DEFAULT_SHARD_SIZE,
            shard_minimizers: DEFAULT_SHARD_MINIMIZERS,
            total_references: 2,
            total_reference_contigs: 2,
            total_mapped_reference_length: 6000,
            total_reference_minimizers: 20,
            total_shard_unique_minimizers: 18,
            build_unix_seconds: 1,
            build_args: vec!["fasterANI".to_string()],
            reference_list_checksum: 42,
            shards: vec![ShardManifestEntry {
                shard_index: 1,
                filename: "database.1.fasketch".to_string(),
                first_reference: 0,
                reference_count: 2,
                reference_contigs: 2,
                mapped_reference_length: 6000,
                reference_minimizers: 20,
                unique_minimizers: 18,
            }],
        }
    }

    #[test]
    fn sharded_sketch_paths_use_prefix_suffixes() {
        let prefix: PathBuf = PathBuf::from("/tmp/database");

        assert_eq!(
            manifest_path(&prefix),
            PathBuf::from("/tmp/database.manifest.json")
        );
        assert_eq!(
            shard_path(&prefix, 2, false),
            PathBuf::from("/tmp/database.2.fasketch")
        );
        assert_eq!(shard_filename(&prefix, 2, false), "database.2.fasketch");
        assert_eq!(
            shard_path(&prefix, 2, true),
            PathBuf::from("/tmp/database.2.fasketch.bgz")
        );
        assert_eq!(shard_filename(&prefix, 2, true), "database.2.fasketch.bgz");
    }

    #[test]
    fn shard_manifest_round_trips_json() -> io::Result<()> {
        let manifest: ShardManifest = sample_shard_manifest();
        let encoded: Vec<u8> = serde_json::to_vec(&manifest).map_err(|err| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("failed to encode test manifest: {err}"),
            )
        })?;
        let decoded: ShardManifest = serde_json::from_slice(&encoded).map_err(|err| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("failed to decode test manifest: {err}"),
            )
        })?;

        assert_eq!(decoded.total_references, manifest.total_references);
        assert_eq!(decoded.shards[0].filename, "database.1.fasketch");

        Ok(())
    }

    #[test]
    fn shard_manifest_rejects_legacy_dust_database() {
        let mut manifest: ShardManifest = sample_shard_manifest();
        manifest.dust_enabled = true;

        let error: Option<String> = shard_manifest_compatibility_error(
            &manifest,
            DEFAULT_KMER_SIZE,
            DEFAULT_WINDOW_SIZE,
            DEFAULT_FRAGMENT_LENGTH,
            DEFAULT_MIN_FRAGMENT_LENGTH,
            DEFAULT_SPLIT_N_RUN,
        );

        assert!(error
            .expect("dust-enabled manifest should be rejected")
            .contains("removed --dust filter"));
    }

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
    fn default_shard_minimizers_scales_with_memory_and_threads() {
        let one_hundred_gib: u64 = 100 * 1024 * 1024 * 1024;
        let shard_minimizers: usize =
            default_shard_minimizers_for_runtime(12, Some(one_hundred_gib));

        assert!(shard_minimizers > DEFAULT_SHARD_MINIMIZERS);
        assert_eq!(
            default_shard_minimizers_for_runtime(12, None),
            DEFAULT_SHARD_MINIMIZERS
        );
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

    #[test]
    fn reference_minimizer_window_estimate_respects_split_n() -> io::Result<()> {
        let path: PathBuf = env::temp_dir().join(format!(
            "fasterani_estimate_windows_{}_{}.fa",
            std::process::id(),
            Instant::now().elapsed().as_nanos()
        ));
        fs::write(&path, b">seq\nAAAAAANNNNAAAAAA\n")?;

        let unsplit_estimate: usize =
            estimate_reference_minimizer_windows(path.to_str().expect("utf8 temp path"), 3, 3, 0)?;
        let split_estimate: usize =
            estimate_reference_minimizer_windows(path.to_str().expect("utf8 temp path"), 3, 3, 4)?;
        fs::remove_file(&path)?;

        assert_eq!(unsplit_estimate, 12);
        assert_eq!(split_estimate, 4);
        assert_eq!(estimate_selected_minimizers_from_windows(12, 3), 6);
        assert_eq!(estimate_selected_minimizers_from_windows(4, 3), 2);

        Ok(())
    }

    #[test]
    fn partitioned_sketch_matches_hash_sketch() -> io::Result<()> {
        let sequence: Vec<u8> = (0..9000)
            .map(|i| b"ACGTGCAATTCG"[i % b"ACGTGCAATTCG".len()])
            .collect();
        let reference_path: PathBuf = env::temp_dir().join(format!(
            "fasterani_partitioned_ref_{}_{}.fa",
            std::process::id(),
            Instant::now().elapsed().as_nanos()
        ));
        let hash_sketch_path: PathBuf = env::temp_dir().join(format!(
            "fasterani_hash_sketch_{}_{}.fasketch",
            std::process::id(),
            Instant::now().elapsed().as_nanos()
        ));
        let partitioned_sketch_path: PathBuf = env::temp_dir().join(format!(
            "fasterani_partitioned_sketch_{}_{}.fasketch",
            std::process::id(),
            Instant::now().elapsed().as_nanos()
        ));
        fs::write(
            &reference_path,
            format!(">ref\n{}\n", String::from_utf8_lossy(&sequence)),
        )?;
        let references: Vec<String> = vec![reference_path.to_string_lossy().into_owned()];

        let hash_stats: SketchBuildStats = ReferenceSketch::collect_and_save_streaming(
            &references,
            DEFAULT_KMER_SIZE,
            DEFAULT_WINDOW_SIZE,
            DEFAULT_FRAGMENT_LENGTH,
            DEFAULT_MIN_FRAGMENT_LENGTH,
            DEFAULT_SPLIT_N_RUN,
            &hash_sketch_path,
            None,
            false,
            1,
            IndexBuildMode::Hash,
            RuntimeOptions::default(),
        )?;
        let partitioned_stats: SketchBuildStats = ReferenceSketch::collect_and_save_streaming(
            &references,
            DEFAULT_KMER_SIZE,
            DEFAULT_WINDOW_SIZE,
            DEFAULT_FRAGMENT_LENGTH,
            DEFAULT_MIN_FRAGMENT_LENGTH,
            DEFAULT_SPLIT_N_RUN,
            &partitioned_sketch_path,
            None,
            false,
            hash_stats.reference_minimizer_count,
            IndexBuildMode::Partitioned,
            RuntimeOptions::default(),
        )?;
        let hash_sketch: ReferenceSketch = ReferenceSketch::load(
            &hash_sketch_path,
            DEFAULT_KMER_SIZE,
            DEFAULT_WINDOW_SIZE,
            DEFAULT_FRAGMENT_LENGTH,
            DEFAULT_MIN_FRAGMENT_LENGTH,
            DEFAULT_SPLIT_N_RUN,
            false,
            None,
            RuntimeOptions::default(),
        )?;
        let partitioned_sketch: ReferenceSketch = ReferenceSketch::load(
            &partitioned_sketch_path,
            DEFAULT_KMER_SIZE,
            DEFAULT_WINDOW_SIZE,
            DEFAULT_FRAGMENT_LENGTH,
            DEFAULT_MIN_FRAGMENT_LENGTH,
            DEFAULT_SPLIT_N_RUN,
            false,
            None,
            RuntimeOptions::default(),
        )?;

        assert_eq!(
            hash_stats.reference_minimizer_count,
            partitioned_stats.reference_minimizer_count
        );
        assert_eq!(
            hash_stats.unique_minimizer_count,
            partitioned_stats.unique_minimizer_count
        );
        assert_eq!(hash_sketch.files.len(), partitioned_sketch.files.len());
        assert_eq!(hash_sketch.contigs.len(), partitioned_sketch.contigs.len());

        for contig_id in 0..hash_sketch.contigs.len() {
            let hash_minimizers: &[ReferenceMinimizer] = hash_sketch
                .contigs
                .minimizers(contig_id)
                .expect("hash contig");
            let partitioned_minimizers: &[ReferenceMinimizer] = partitioned_sketch
                .contigs
                .minimizers(contig_id)
                .expect("partitioned contig");
            assert_eq!(hash_minimizers, partitioned_minimizers);

            for minimizer in hash_minimizers {
                assert_eq!(
                    hash_sketch.index.get(&minimizer.hash),
                    partitioned_sketch.index.get(&minimizer.hash)
                );
            }
        }

        fs::remove_file(reference_path)?;
        fs::remove_file(contig_sidecar_path(&hash_sketch_path))?;
        fs::remove_file(contig_sidecar_path(&partitioned_sketch_path))?;
        fs::remove_file(hash_sketch_path)?;
        fs::remove_file(partitioned_sketch_path)?;

        Ok(())
    }

    #[test]
    fn minimizers_skip_kmers_with_ambiguous_bases() {
        let sequence: &[u8] =
            b"AGCTTAGGCTAACCGTATGCCGATTAACGNNNNNNNNNNGCTAGTCCATGATCGTACCGTTAAGGCTA";
        let kmer_size: usize = 10usize;
        let window_size: usize = 16usize;

        let minimizers: Vec<(MinimizerKey, u32)> =
            canonical_minimizers_with_positions(sequence, kmer_size, window_size);

        assert!(!minimizers.is_empty());
        for (_hash, position) in minimizers {
            let position: usize = position as usize;
            let kmer: &[u8] = &sequence[position..position + kmer_size];
            assert!(!kmer.iter().any(|base| matches!(base, b'N' | b'n')));
        }
    }

    #[test]
    fn all_ambiguous_sequence_has_no_minimizers() {
        let sequence: &[u8] = b"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
        let kmer_size: usize = 10usize;
        let window_size: usize = 16usize;

        let minimizers: Vec<(MinimizerKey, u32)> =
            canonical_minimizers_with_positions(sequence, kmer_size, window_size);

        assert!(minimizers.is_empty());
    }

    #[test]
    fn expected_minimizer_windows_match_fragment_geometry() {
        assert_eq!(expected_minimizer_window_count(3000, 16, 24), 2962);
        assert_eq!(expected_minimizer_window_count(38, 16, 24), 0);
        assert_eq!(expected_minimizer_window_count(39, 16, 24), 1);
    }

    #[test]
    fn default_fragment_settings_are_fastani_compatible() {
        assert_eq!(DEFAULT_FRAGMENT_STRIDE, DEFAULT_FRAGMENT_LENGTH);
        assert_eq!(DEFAULT_MIN_FRAGMENT_LENGTH, DEFAULT_FRAGMENT_LENGTH);
        assert!(fastani_compatible_fragment_mode(
            DEFAULT_FRAGMENT_LENGTH,
            DEFAULT_FRAGMENT_STRIDE,
            DEFAULT_MIN_FRAGMENT_LENGTH
        ));
    }

    #[test]
    fn clean_sequence_has_no_missing_minimizer_windows() {
        let sequence: Vec<u8> = repeated_acgt(120);

        assert_eq!(
            usable_minimizer_window_count(&sequence, 5, 5),
            expected_minimizer_window_count(sequence.len(), 5, 5)
        );
    }

    #[test]
    fn all_ambiguous_observation_keeps_expected_window_count() {
        let sequence: &[u8] = b"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
        let observation: MinimizerObservation = canonical_minimizer_observation(sequence, 10, 16);

        assert_eq!(usable_minimizer_window_count(sequence, 10, 16), 0);
        assert_eq!(expected_minimizer_window_count(sequence.len(), 10, 16), 16);
        assert!(observation.minimizers_with_positions.is_empty());
    }

    #[test]
    fn minmer_seeds_use_smallest_hashes() {
        let minimizers: Vec<MinimizerKey> = vec![3, 5, 8, 13, 21];

        let seeds: Vec<MinimizerKey> = select_seed_minimizers(&minimizers, Some(3));

        assert_eq!(seeds, vec![3, 5, 8]);
    }

    #[test]
    fn sketch_reference_name_keeps_only_basename() {
        assert_eq!(
            sketch_reference_name("/tmp/reference/GCF_000146045.2_R64_genomic.fna"),
            "GCF_000146045.2_R64_genomic.fna"
        );
        assert_eq!(sketch_reference_name("relative.fa"), "relative.fa");
    }

    #[test]
    fn loaded_sketch_contig_minimizers_match_owned_minimizers() -> io::Result<()> {
        let contigs: Vec<ReferenceContig> = vec![
            ReferenceContig {
                file_id: 0,
                minimizers: vec![
                    ReferenceMinimizer {
                        hash: 11,
                        position: 3,
                    },
                    ReferenceMinimizer {
                        hash: 17,
                        position: 9,
                    },
                ],
            },
            ReferenceContig {
                file_id: 0,
                minimizers: vec![ReferenceMinimizer {
                    hash: 23,
                    position: 4,
                }],
            },
        ];
        let mut index: ReferenceHitMap = ReferenceHitMap::default();
        for (contig_id, contig) in contigs.iter().enumerate() {
            for minimizer in &contig.minimizers {
                index.entry(minimizer.hash).or_default().push(SeedHit {
                    reference_contig_id: contig_id as u32,
                    position: minimizer.position,
                });
            }
        }

        let sketch: ReferenceSketch = ReferenceSketch {
            files: vec![ReferenceFile {
                path: "/tmp/ref.fa".to_string(),
                mapped_length: 3000,
            }],
            contigs: ReferenceContigs::Owned(contigs.clone()),
            contig_names: Some(vec![
                ReferenceContigName {
                    file_id: 0,
                    name: "ref_contig_a".to_string(),
                    segment_start: 0,
                    segment_end: 100,
                },
                ReferenceContigName {
                    file_id: 0,
                    name: "ref_contig_b".to_string(),
                    segment_start: 200,
                    segment_end: 300,
                },
            ]),
            index: ReferenceIndex::Hash(index),
        };
        let path: PathBuf = env::temp_dir().join(format!(
            "fasterani_mmap_contigs_{}_{}.fasketch",
            std::process::id(),
            Instant::now().elapsed().as_nanos()
        ));

        sketch.save(
            &path,
            DEFAULT_KMER_SIZE,
            DEFAULT_WINDOW_SIZE,
            DEFAULT_FRAGMENT_LENGTH,
            DEFAULT_MIN_FRAGMENT_LENGTH,
            0,
            None,
            RuntimeOptions::default(),
        )?;
        let loaded: ReferenceSketch = ReferenceSketch::load(
            &path,
            DEFAULT_KMER_SIZE,
            DEFAULT_WINDOW_SIZE,
            DEFAULT_FRAGMENT_LENGTH,
            DEFAULT_MIN_FRAGMENT_LENGTH,
            0,
            true,
            None,
            RuntimeOptions::default(),
        )?;
        let loaded_contig_names: &[ReferenceContigName] =
            loaded.contig_names.as_deref().expect("loaded sidecar");
        assert_eq!(loaded_contig_names[0].name, "ref_contig_a");
        assert_eq!(loaded_contig_names[1].segment_start, 200);
        fs::remove_file(&path)?;
        fs::remove_file(contig_sidecar_path(&path))?;

        assert_eq!(loaded.files[0].path, "ref.fa");
        assert_eq!(loaded.contigs.len(), contigs.len());
        for (contig_id, contig) in contigs.iter().enumerate() {
            assert_eq!(loaded.contigs.file_id(contig_id), Some(contig.file_id));
            assert_eq!(
                loaded.contigs.minimizers(contig_id).expect("loaded contig"),
                contig.minimizers.as_slice()
            );
        }

        Ok(())
    }

    #[test]
    fn default_fragment_ranges_match_fastani_chunks() {
        let ranges: Vec<std::ops::Range<usize>> = query_fragment_ranges(7500, 3000, 3000, 3000);

        assert_eq!(ranges, vec![0..3000, 3000..6000]);
    }

    #[test]
    fn default_fragment_ranges_discard_short_terminal_remainder() {
        let ranges: Vec<std::ops::Range<usize>> = query_fragment_ranges(8999, 3000, 3000, 3000);

        assert_eq!(ranges, vec![0..3000, 3000..6000]);
    }

    #[test]
    fn adaptive_fragment_ranges_overlap_and_cover_tail() {
        let ranges: Vec<std::ops::Range<usize>> = query_fragment_ranges(7500, 3000, 1000, 1000);

        assert_eq!(
            ranges,
            vec![
                0..3000,
                1000..4000,
                2000..5000,
                3000..6000,
                4000..7000,
                5000..7500
            ]
        );
    }

    #[test]
    fn adaptive_fragment_ranges_keep_short_usable_contigs() {
        let ranges: Vec<std::ops::Range<usize>> = query_fragment_ranges(1500, 3000, 1000, 1000);

        assert_eq!(ranges, vec![0..1500]);
    }

    #[test]
    fn adaptive_fragment_ranges_keep_terminal_remainder() {
        let ranges: Vec<std::ops::Range<usize>> = query_fragment_ranges(7500, 3000, 3000, 1000);

        assert_eq!(ranges, vec![0..3000, 3000..6000, 6000..7500]);
    }

    #[test]
    fn mapped_length_keeps_terminal_remainder_when_allowed() {
        assert_eq!(
            mapped_length_from_fragment_ranges(7500, DEFAULT_FRAGMENT_LENGTH, 3000),
            6000
        );
        assert_eq!(
            mapped_length_from_fragment_ranges(7500, DEFAULT_FRAGMENT_LENGTH, 1000),
            7500
        );
        assert_eq!(
            mapped_length_from_fragment_ranges(3500, DEFAULT_FRAGMENT_LENGTH, 1000),
            3000
        );
        assert_eq!(
            mapped_length_from_fragment_ranges(4000, DEFAULT_FRAGMENT_LENGTH, 1000),
            4000
        );
    }

    #[test]
    fn split_n_zero_keeps_whole_sequence() {
        let ranges: Vec<std::ops::Range<usize>> = split_sequence_ranges(b"ACGTNNNNACGT", 0);

        assert_eq!(ranges, vec![0..12]);
    }

    #[test]
    fn split_n_keeps_short_ambiguous_runs() {
        let ranges: Vec<std::ops::Range<usize>> = split_sequence_ranges(b"ACGTNNACGT", 3);

        assert_eq!(ranges, vec![0..10]);
    }

    #[test]
    fn split_n_breaks_long_ambiguous_runs() {
        let ranges: Vec<std::ops::Range<usize>> = split_sequence_ranges(b"NNNACGTNNNNACGTNNN", 3);

        assert_eq!(ranges, vec![3..7, 11..15]);
    }
}

fn fastani_mash_distance(jaccard: f64, kmer_size: usize) -> f64 {
    if jaccard == 0.0 {
        return 1.0;
    }

    if jaccard == 1.0 {
        return 0.0;
    }

    (-1.0 / kmer_size as f64) * ((2.0 * jaccard) / (1.0 + jaccard)).ln()
}

fn mash_distance_to_jaccard(distance: f64, kmer_size: usize) -> f64 {
    1.0 / (2.0 * (kmer_size as f64 * distance).exp() - 1.0)
}

fn binomial_survival_at_least(x: usize, probability: f64, trials: usize) -> f64 {
    if x == 0 {
        return 1.0;
    }

    if x > trials || probability <= 0.0 {
        return 0.0;
    }

    if probability >= 1.0 {
        return 1.0;
    }

    let q = 1.0 - probability;
    let mut pmf = q.powi(trials as i32);
    let mut cdf_below = pmf;

    for i in 0..(x - 1) {
        pmf *= (trials - i) as f64 / (i + 1) as f64 * probability / q;
        cdf_below += pmf;
    }

    (1.0 - cdf_below).clamp(0.0, 1.0)
}

fn mash_distance_lower_bound(
    distance: f64,
    sketch_size: usize,
    kmer_size: usize,
    confidence_interval: f64,
) -> f64 {
    let q2: f64 = (1.0 - confidence_interval) / 2.0;
    let jaccard: f64 = mash_distance_to_jaccard(distance, kmer_size);
    let mut x: usize = ((sketch_size as f64 * jaccard).ceil() as usize).max(1);

    while x <= sketch_size {
        let cdf_complement = binomial_survival_at_least(x, jaccard, sketch_size);

        if cdf_complement < q2 {
            x = x.saturating_sub(1);
            break;
        }

        x += 1;
    }

    fastani_mash_distance(x as f64 / sketch_size as f64, kmer_size)
}

fn estimate_minimum_shared_minimizers(
    sketch_size: usize,
    kmer_size: usize,
    percent_identity: f64,
) -> usize {
    let mash_distance: f64 = 1.0 - percent_identity / 100.0;
    let jaccard: f64 = mash_distance_to_jaccard(mash_distance, kmer_size);

    (sketch_size as f64 * jaccard).ceil() as usize
}

/// Estimate the relaxed minimum shared minimizers needed to keep a candidate alive.
fn estimate_relaxed_minimum_shared_minimizers(
    sketch_size: usize,
    kmer_size: usize,
    percent_identity: f64,
    mash_confidence: f64,
) -> usize {
    let strict_minimum: usize =
        estimate_minimum_shared_minimizers(sketch_size, kmer_size, percent_identity);
    let mut relaxed_minimum: usize = strict_minimum;

    for i in (0..=strict_minimum).rev() {
        let jaccard = i as f64 / sketch_size as f64;
        let distance = fastani_mash_distance(jaccard, kmer_size);
        let lower_distance =
            mash_distance_lower_bound(distance, sketch_size, kmer_size, mash_confidence);
        let upper_identity = 100.0 * (1.0 - lower_distance);

        if upper_identity >= percent_identity {
            relaxed_minimum = i;
        } else {
            break;
        }
    }

    relaxed_minimum
}

fn lower_bound_minimizer_position(minimizers: &[ReferenceMinimizer], position: u32) -> usize {
    minimizers.partition_point(|minimizer| minimizer.position < position)
}

fn align_up(value: usize, alignment: usize) -> usize {
    debug_assert!(alignment.is_power_of_two());
    (value + alignment - 1) & !(alignment - 1)
}

fn checked_section_end(offset: usize, count: usize, item_size: usize) -> io::Result<usize> {
    let byte_len = count.checked_mul(item_size).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            "reference sketch section size overflow",
        )
    })?;
    offset.checked_add(byte_len).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            "reference sketch section offset overflow",
        )
    })
}

fn sketch_reference_name(path: &str) -> String {
    Path::new(path)
        .file_name()
        .map(|name| name.to_string_lossy().into_owned())
        .unwrap_or_else(|| path.to_string())
}

fn append_path_suffix(path: &Path, suffix: &str) -> PathBuf {
    let mut value = path.as_os_str().to_os_string();
    value.push(suffix);
    PathBuf::from(value)
}

fn path_has_extension(path: &Path, expected: &str) -> bool {
    path.extension()
        .and_then(|extension| extension.to_str())
        .is_some_and(|extension| extension.eq_ignore_ascii_case(expected))
}

fn is_gzip_path(path: &Path) -> bool {
    path_has_extension(path, "gz") || path_has_extension(path, "bgz")
}

fn gzp_error_to_io(error: gzp::GzpError) -> io::Error {
    io::Error::new(
        io::ErrorKind::Other,
        format!("failed to finish BGZF compression: {error}"),
    )
}

fn open_fasta_reader(path: &str) -> io::Result<fasta::io::Reader<Box<dyn io::BufRead>>> {
    let path_ref: &Path = Path::new(path);
    let file: fs::File = fs::File::open(path_ref)?;
    let reader: Box<dyn io::BufRead> = if is_gzip_path(path_ref) {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    fasta::io::reader::Builder.build_from_reader(reader)
}

fn read_text_maybe_gzip(path: &Path) -> io::Result<String> {
    let file: fs::File = fs::File::open(path)?;
    let mut contents: String = String::new();

    if is_gzip_path(path) {
        let mut reader: BufReader<MultiGzDecoder<fs::File>> =
            BufReader::new(MultiGzDecoder::new(file));
        reader.read_to_string(&mut contents)?;
    } else {
        let mut reader: BufReader<fs::File> = BufReader::new(file);
        reader.read_to_string(&mut contents)?;
    }

    Ok(contents)
}

fn compress_file_to_bgzf(source: &Path, destination: &Path, threads: usize) -> io::Result<()> {
    if let Some(parent) = destination
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        fs::create_dir_all(parent)?;
    }

    let input: fs::File = fs::File::open(source)?;
    let output: fs::File = fs::File::create(destination)?;
    let mut reader: BufReader<fs::File> = BufReader::new(input);
    let mut writer = ZBuilder::<Bgzf, fs::File>::new()
        .num_threads(threads.max(1))
        .from_writer(output);
    io::copy(&mut reader, &mut writer)?;
    writer.finish().map_err(gzp_error_to_io)?;

    Ok(())
}

fn decompress_to_scratch(
    source: &Path,
    tmp_dir: Option<&Path>,
    purpose: &str,
) -> io::Result<ScratchFile> {
    let (scratch, scratch_file): (ScratchFile, fs::File) = ScratchFile::create(tmp_dir, purpose)?;
    let input: fs::File = fs::File::open(source)?;
    let mut reader: BufReader<MultiGzDecoder<fs::File>> =
        BufReader::new(MultiGzDecoder::new(input));
    let mut writer: BufWriter<fs::File> = BufWriter::new(scratch_file);
    io::copy(&mut reader, &mut writer)?;
    writer.flush()?;
    drop(writer);

    Ok(scratch)
}

struct SketchOutput {
    final_path: PathBuf,
    write_path: PathBuf,
    scratch: Option<ScratchFile>,
    writer: Option<BufWriter<fs::File>>,
    bgzip: bool,
    threads: usize,
}

impl SketchOutput {
    fn create(
        final_path: &Path,
        tmp_dir: Option<&Path>,
        bgzip: bool,
        threads: usize,
    ) -> io::Result<Self> {
        if let Some(parent) = final_path
            .parent()
            .filter(|parent| !parent.as_os_str().is_empty())
        {
            fs::create_dir_all(parent)?;
        }

        if bgzip {
            let (scratch, file): (ScratchFile, fs::File) =
                ScratchFile::create(tmp_dir, "uncompressed-sketch")?;
            let write_path: PathBuf = scratch.path.clone();
            Ok(Self {
                final_path: final_path.to_path_buf(),
                write_path,
                scratch: Some(scratch),
                writer: Some(BufWriter::new(file)),
                bgzip,
                threads,
            })
        } else {
            let file: fs::File = fs::File::create(final_path)?;
            Ok(Self {
                final_path: final_path.to_path_buf(),
                write_path: final_path.to_path_buf(),
                scratch: None,
                writer: Some(BufWriter::new(file)),
                bgzip,
                threads,
            })
        }
    }

    fn writer_mut(&mut self) -> io::Result<&mut BufWriter<fs::File>> {
        self.writer.as_mut().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::Other,
                "reference sketch output writer is already closed",
            )
        })
    }

    fn finish(mut self) -> io::Result<u64> {
        if let Some(mut writer) = self.writer.take() {
            writer.flush()?;
        }

        if self.bgzip {
            compress_file_to_bgzf(&self.write_path, &self.final_path, self.threads)?;
        }

        let output_len: u64 = fs::metadata(&self.final_path)?.len();
        drop(self.scratch.take());
        Ok(output_len)
    }
}

fn contig_sidecar_path(sketch_path: &Path) -> PathBuf {
    let sketch_path_string: String = sketch_path.to_string_lossy().into_owned();
    if let Some(uncompressed_name) = sketch_path_string.strip_suffix(".bgz") {
        return append_path_suffix(Path::new(uncompressed_name), ".contigs.tsv.bgz");
    }

    append_path_suffix(sketch_path, ".contigs.tsv")
}

fn manifest_path(prefix: &Path) -> PathBuf {
    append_path_suffix(prefix, ".manifest.json")
}

fn shard_path(prefix: &Path, shard_index: usize, bgzip: bool) -> PathBuf {
    if bgzip {
        append_path_suffix(prefix, &format!(".{shard_index}.fasketch.bgz"))
    } else {
        append_path_suffix(prefix, &format!(".{shard_index}.fasketch"))
    }
}

fn shard_filename(prefix: &Path, shard_index: usize, bgzip: bool) -> String {
    shard_path(prefix, shard_index, bgzip)
        .file_name()
        .map(|name| name.to_string_lossy().into_owned())
        .unwrap_or_else(|| {
            shard_path(prefix, shard_index, bgzip)
                .to_string_lossy()
                .into_owned()
        })
}

fn shard_entry_path(prefix: &Path, entry: &ShardManifestEntry) -> PathBuf {
    let filename_path: PathBuf = PathBuf::from(&entry.filename);
    if filename_path.is_absolute() {
        return filename_path;
    }

    prefix
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .map(|parent| parent.join(&entry.filename))
        .unwrap_or(filename_path)
}

fn legacy_sketch_path(prefix: &Path) -> Option<PathBuf> {
    if prefix.exists() {
        return Some(prefix.to_path_buf());
    }

    if prefix.extension().is_none() {
        for extension in ["fasketch.bgz", "fasketch"] {
            let path: PathBuf = append_path_suffix(prefix, &format!(".{extension}"));
            if path.exists() {
                return Some(path);
            }
        }
    }

    None
}

fn reference_list_checksum(reference_paths: &[String]) -> u64 {
    const FNV_OFFSET: u64 = 0xcbf29ce484222325;
    const FNV_PRIME: u64 = 0x100000001b3;

    let mut hash: u64 = FNV_OFFSET;
    for path in reference_paths {
        for byte in path.as_bytes().iter().copied().chain(std::iter::once(0)) {
            hash ^= u64::from(byte);
            hash = hash.wrapping_mul(FNV_PRIME);
        }
    }

    hash
}

fn tsv_field(value: &str) -> String {
    value
        .chars()
        .map(|ch| match ch {
            '\t' | '\n' | '\r' => ' ',
            _ => ch,
        })
        .collect()
}

fn write_contig_name_sidecar(
    sketch_path: &Path,
    files: &[ReferenceFile],
    contig_names: &[ReferenceContigName],
    threads: usize,
) -> io::Result<()> {
    let path: PathBuf = contig_sidecar_path(sketch_path);
    if let Some(parent) = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        fs::create_dir_all(parent)?;
    }

    let write_records = |writer: &mut dyn Write| -> io::Result<()> {
        writeln!(
            writer,
            "contig_id\treference_file_id\treference_file\treference_contig\tsegment_start\tsegment_end"
        )?;

        for (contig_id, contig) in contig_names.iter().enumerate() {
            let reference_file: &str = files
                .get(contig.file_id)
                .map(|file| file.path.as_str())
                .unwrap_or("unknown");
            writeln!(
                writer,
                "{contig_id}\t{}\t{}\t{}\t{}\t{}",
                contig.file_id,
                tsv_field(reference_file),
                tsv_field(&contig.name),
                contig.segment_start,
                contig.segment_end
            )?;
        }

        Ok(())
    };

    if is_gzip_path(&path) {
        let file: fs::File = fs::File::create(path)?;
        let mut writer = ZBuilder::<Bgzf, fs::File>::new()
            .num_threads(threads.max(1))
            .from_writer(file);
        write_records(&mut *writer)?;
        writer.finish().map_err(gzp_error_to_io)?;
    } else {
        let file: fs::File = fs::File::create(path)?;
        let mut writer: BufWriter<fs::File> = BufWriter::new(file);
        write_records(&mut writer)?;
        writer.flush()?;
    }

    Ok(())
}

fn load_contig_sidecar_contents(path: &Path) -> io::Result<String> {
    read_text_maybe_gzip(path).map_err(|err| {
        io::Error::new(
            err.kind(),
            format!(
                "failed to read contig sidecar {}; rebuild the sketch database to create it: {err}",
                path.display()
            ),
        )
    })
}

fn uncompressed_contig_sidecar_path_for_bgzip_sketch(sketch_path: &Path) -> Option<PathBuf> {
    let sketch_path_string: String = sketch_path.to_string_lossy().into_owned();
    sketch_path_string
        .strip_suffix(".bgz")
        .map(|uncompressed_name| append_path_suffix(Path::new(uncompressed_name), ".contigs.tsv"))
}

fn read_contig_sidecar_text(sketch_path: &Path) -> io::Result<(PathBuf, String)> {
    let path: PathBuf = contig_sidecar_path(sketch_path);
    match load_contig_sidecar_contents(&path) {
        Ok(contents) => Ok((path, contents)),
        Err(error) if error.kind() == io::ErrorKind::NotFound => {
            if let Some(fallback_path) =
                uncompressed_contig_sidecar_path_for_bgzip_sketch(sketch_path)
            {
                let contents: String = load_contig_sidecar_contents(&fallback_path)?;
                Ok((fallback_path, contents))
            } else {
                Err(error)
            }
        }
        Err(error) => Err(error),
    }
}

fn parse_usize_field(value: &str, field: &str, path: &Path) -> io::Result<usize> {
    value.parse::<usize>().map_err(|err| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "failed to parse {field} in contig sidecar {}: {err}",
                path.display()
            ),
        )
    })
}

fn parse_u32_field(value: &str, field: &str, path: &Path) -> io::Result<u32> {
    value.parse::<u32>().map_err(|err| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "failed to parse {field} in contig sidecar {}: {err}",
                path.display()
            ),
        )
    })
}

fn load_contig_name_sidecar(
    sketch_path: &Path,
    expected_contigs: usize,
) -> io::Result<Vec<ReferenceContigName>> {
    let (path, contents): (PathBuf, String) = read_contig_sidecar_text(sketch_path)?;
    let mut lines = contents.lines();
    let header: &str = lines.next().ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("contig sidecar {} is empty", path.display()),
        )
    })?;
    if header
        != "contig_id\treference_file_id\treference_file\treference_contig\tsegment_start\tsegment_end"
    {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("contig sidecar {} has an unexpected header", path.display()),
        ));
    }

    let mut contigs: Vec<ReferenceContigName> = Vec::with_capacity(expected_contigs);
    for (line_index, line) in lines.enumerate() {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 6 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "contig sidecar {} line {} has {} fields, expected 6",
                    path.display(),
                    line_index + 2,
                    fields.len()
                ),
            ));
        }

        let contig_id: usize = parse_usize_field(fields[0], "contig_id", &path)?;
        if contig_id != contigs.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "contig sidecar {} expected contig_id {}, found {contig_id}",
                    path.display(),
                    contigs.len()
                ),
            ));
        }

        contigs.push(ReferenceContigName {
            file_id: parse_usize_field(fields[1], "reference_file_id", &path)?,
            name: fields[3].to_string(),
            segment_start: parse_u32_field(fields[4], "segment_start", &path)?,
            segment_end: parse_u32_field(fields[5], "segment_end", &path)?,
        });
    }

    if contigs.len() != expected_contigs {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "contig sidecar {} has {} contigs, expected {expected_contigs}",
                path.display(),
                contigs.len()
            ),
        ));
    }

    Ok(contigs)
}

fn unix_timestamp_seconds() -> io::Result<u64> {
    Ok(SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_err(|err| {
            io::Error::new(
                io::ErrorKind::Other,
                format!("system clock is before UNIX epoch: {err}"),
            )
        })?
        .as_secs())
}

fn validate_shard_size(shard_size: usize) -> io::Result<()> {
    if shard_size == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--shard-size must be at least 1",
        ));
    }

    Ok(())
}

fn default_fragment_length() -> u32 {
    DEFAULT_FRAGMENT_LENGTH
}

fn default_shard_minimizers() -> usize {
    DEFAULT_SHARD_MINIMIZERS
}

fn default_shard_minimizers_for_runtime(threads: usize, max_memory_bytes: Option<u64>) -> usize {
    let Some(max_memory_bytes) = max_memory_bytes else {
        return DEFAULT_SHARD_MINIMIZERS;
    };

    let Ok(max_memory_bytes) = usize::try_from(max_memory_bytes) else {
        return DEFAULT_SHARD_MINIMIZERS;
    };

    let active_jobs: usize = threads.max(1);
    let per_job_bytes: usize = max_memory_bytes / active_jobs;
    let memory_sized_minimizers: usize =
        per_job_bytes / ESTIMATED_PARTITIONED_SHARD_BYTES_PER_MINIMIZER;

    memory_sized_minimizers.max(DEFAULT_SHARD_MINIMIZERS)
}

fn validate_shard_minimizers(shard_minimizers: usize) -> io::Result<()> {
    if shard_minimizers == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--shard-minimizers must be at least 1",
        ));
    }

    Ok(())
}

fn validate_kmer_size(kmer_size: usize) -> io::Result<()> {
    if !(1..=16).contains(&kmer_size) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--kmer-size must be between 1 and 16 for canonical-2bit-u32 keys",
        ));
    }

    Ok(())
}

fn validate_window_size(window_size: usize) -> io::Result<()> {
    if window_size == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--window-size must be at least 1",
        ));
    }

    Ok(())
}

fn validate_fragment_length(fragment_length: u32) -> io::Result<()> {
    if fragment_length == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--fragment-length must be at least 1",
        ));
    }

    Ok(())
}

fn validate_min_identity(min_identity: f64) -> io::Result<()> {
    if !min_identity.is_finite() || !(0.0..=100.0).contains(&min_identity) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--min-identity must be a finite value between 0 and 100",
        ));
    }

    Ok(())
}

fn validate_mash_confidence(mash_confidence: f64) -> io::Result<()> {
    if !mash_confidence.is_finite() || !(0.0..1.0).contains(&mash_confidence) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--mash-confidence must be a finite value in [0, 1)",
        ));
    }

    Ok(())
}

fn estimate_partitioned_shard_memory_bytes(estimated_minimizers: usize) -> usize {
    estimated_minimizers
        .saturating_mul(ESTIMATED_PARTITIONED_SHARD_BYTES_PER_MINIMIZER)
        .max(DEFAULT_PARTITION_TARGET_BYTES)
}

fn database_build_parallelism(
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

fn partition_build_plan(
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

fn effective_index_build_mode(
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

fn partition_id_for_key(key: MinimizerKey, partition_count: usize) -> usize {
    debug_assert!(partition_count.is_power_of_two());
    let partition_bits: u32 = partition_count.trailing_zeros();
    if partition_bits == 0 {
        return 0;
    }

    (key >> (MinimizerKey::BITS - partition_bits)) as usize
}

fn estimate_reference_minimizer_windows(
    reference_path: &str,
    kmer_size: usize,
    window_size: usize,
    split_n_run: usize,
) -> io::Result<usize> {
    let mut reader: fasta::io::Reader<Box<dyn io::BufRead>> = open_fasta_reader(reference_path)?;
    let mut minimizer_window_count: usize = 0usize;

    for result in reader.records() {
        let record: fasta::Record = result.map_err(|err| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("failed to read FASTA record from reference {reference_path}: {err}"),
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

fn estimate_selected_minimizers_from_windows(window_count: usize, window_size: usize) -> usize {
    let denominator: usize = window_size.saturating_add(1);
    if denominator == 0 {
        return window_count;
    }

    window_count
        .saturating_mul(2)
        .saturating_add(denominator - 1)
        / denominator
}

fn plan_shards_from_minimizer_counts(
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

fn plan_shards_by_minimizers(
    reference_paths: &[String],
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
    let planner_threads: usize = threads.max(1).min(reference_paths.len().max(1));

    if runtime_options.progress_enabled {
        emit_progress(
            "shard_plan",
            &format!(
                "event=start\testimator=window-count\treferences={}\tshard_size={shard_size}\tshard_minimizers={shard_minimizers}\tthreads={threads}\tplanner_parallelism={planner_threads}",
                reference_paths.len(),
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
    let count_results: Vec<io::Result<(usize, usize)>> = pool.install(|| {
        reference_paths
            .par_iter()
            .enumerate()
            .map(|(reference_index, reference_path)| {
                let reference_minimizer_windows: usize = estimate_reference_minimizer_windows(
                    reference_path,
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
                    && (references_done % REFERENCE_PROGRESS_INTERVAL == 0
                        || references_done == reference_paths.len())
                {
                    emit_progress(
                        "shard_plan",
                        &format!(
                            "event=references\treferences_done={references_done}\treferences_total={}\tplanner_parallelism={planner_threads}",
                            reference_paths.len(),
                        ),
                        plan_start,
                    );
                }
                check_memory_limit("during shard planning", runtime_options)?;

                Ok((reference_index, reference_minimizers))
            })
            .collect()
    });

    let mut minimizer_counts: Vec<usize> = vec![0usize; reference_paths.len()];
    for result in count_results {
        let (reference_index, reference_minimizers): (usize, usize) = result?;
        minimizer_counts[reference_index] = reference_minimizers;
    }

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

fn shard_manifest_compatibility_error(
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

fn slice_as_bytes<T>(slice: &[T]) -> &[u8] {
    unsafe { std::slice::from_raw_parts(slice.as_ptr() as *const u8, std::mem::size_of_val(slice)) }
}

fn slice_as_bytes_mut<T>(slice: &mut [T]) -> &mut [u8] {
    unsafe {
        std::slice::from_raw_parts_mut(slice.as_mut_ptr() as *mut u8, std::mem::size_of_val(slice))
    }
}

fn write_padding(writer: &mut impl Write, len: usize) -> io::Result<()> {
    const ZEROES: [u8; 8] = [0; 8];
    writer.write_all(&ZEROES[..len])
}

struct ScratchFile {
    path: PathBuf,
}

impl ScratchFile {
    fn create(tmp_dir: Option<&Path>, purpose: &str) -> io::Result<(Self, fs::File)> {
        let directory: PathBuf = tmp_dir.map(Path::to_path_buf).unwrap_or_else(env::temp_dir);
        fs::create_dir_all(&directory)?;

        for attempt in 0..100u32 {
            let timestamp: u128 = SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .map_err(|err| {
                    io::Error::new(
                        io::ErrorKind::Other,
                        format!("system clock is before UNIX epoch: {err}"),
                    )
                })?
                .as_nanos();
            let path: PathBuf = directory.join(format!(
                "fasterani-{}-{timestamp}-{purpose}-{attempt}.tmp",
                std::process::id()
            ));
            match fs::OpenOptions::new()
                .create_new(true)
                .read(true)
                .write(true)
                .open(&path)
            {
                Ok(file) => return Ok((Self { path }, file)),
                Err(error) if error.kind() == io::ErrorKind::AlreadyExists => continue,
                Err(error) => return Err(error),
            }
        }

        Err(io::Error::new(
            io::ErrorKind::AlreadyExists,
            "failed to create a unique fasterANI scratch file",
        ))
    }
}

impl Drop for ScratchFile {
    fn drop(&mut self) {
        let _ = fs::remove_file(&self.path);
    }
}

struct PartitionWriters {
    files: Vec<ScratchFile>,
    writers: Vec<BufWriter<fs::File>>,
    buffers: Vec<Vec<PartitionHitRecord>>,
    buffer_record_limit: usize,
}

impl PartitionWriters {
    fn new(
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

    fn partition_count(&self) -> usize {
        self.files.len()
    }

    fn path(&self, partition_index: usize) -> &Path {
        &self.files[partition_index].path
    }

    fn push(&mut self, record: PartitionHitRecord) -> io::Result<()> {
        let partition_index: usize = partition_id_for_key(record.key, self.partition_count());
        let buffer: &mut Vec<PartitionHitRecord> = &mut self.buffers[partition_index];
        buffer.push(record);

        if buffer.len() >= self.buffer_record_limit {
            self.flush_partition(partition_index)?;
        }

        Ok(())
    }

    fn flush_partition(&mut self, partition_index: usize) -> io::Result<()> {
        let buffer: &mut Vec<PartitionHitRecord> = &mut self.buffers[partition_index];
        if buffer.is_empty() {
            return Ok(());
        }

        self.writers[partition_index].write_all(slice_as_bytes(buffer))?;
        buffer.clear();

        Ok(())
    }

    fn flush_all(&mut self) -> io::Result<()> {
        for partition_index in 0..self.partition_count() {
            self.flush_partition(partition_index)?;
            self.writers[partition_index].flush()?;
        }

        Ok(())
    }
}

impl ReferenceSketch {
    #[cfg(debug_assertions)]
    fn memory_estimate(&self) -> ReferenceMemoryEstimate {
        let reference_minimizers: usize = self.contigs.total_minimizers();
        let reference_minimizer_vec_bytes: usize = self.contigs.owned_minimizer_capacity_bytes();

        match &self.index {
            ReferenceIndex::Hash(index) => {
                let seed_hits: usize = index.values().map(Vec::len).sum();
                let seed_hit_vec_bytes: usize = index
                    .values()
                    .map(|hits| hits.capacity() * size_of::<SeedHit>())
                    .sum();
                let hash_index_rough_bytes: usize =
                    index.capacity() * (size_of::<MinimizerKey>() + size_of::<Vec<SeedHit>>() + 8);

                ReferenceMemoryEstimate {
                    reference_minimizers,
                    reference_minimizer_vec_bytes,
                    unique_index_keys: index.len(),
                    seed_hits,
                    seed_hit_vec_bytes,
                    hash_index_rough_bytes,
                    ..ReferenceMemoryEstimate::default()
                }
            }
            ReferenceIndex::Mphf(index) => ReferenceMemoryEstimate {
                reference_minimizers,
                reference_minimizer_vec_bytes,
                unique_index_keys: index.key_count,
                seed_hits: index.hit_count,
                mmap_file_bytes: index.mmap.len,
                mmap_slot_key_bytes: index.key_count * size_of::<MinimizerKey>(),
                mmap_hit_offset_bytes: index.key_count * size_of::<u64>(),
                mmap_hit_count_bytes: index.key_count * size_of::<u32>(),
                mmap_hit_payload_bytes: index.hit_count * size_of::<SeedHit>(),
                mmap_contig_record_bytes: self.contigs.len() * size_of::<ContigRecord>(),
                mmap_reference_minimizer_bytes: reference_minimizers
                    * size_of::<ReferenceMinimizer>(),
                ..ReferenceMemoryEstimate::default()
            },
        }
    }

    /// Build a reference sketch from all provided reference FASTA files.
    fn collect(
        reference_paths: &[String],
        kmer_size: usize,
        window_size: usize,
        fragment_length: u32,
        min_fragment_length: u32,
        split_n_run: usize,
        runtime_options: RuntimeOptions,
    ) -> io::Result<Self> {
        let mut files: Vec<ReferenceFile> = Vec::new();
        let mut contigs: Vec<ReferenceContig> = Vec::new();
        let mut index: ReferenceHitMap = ReferenceHitMap::default();
        let build_start: Instant = Instant::now();
        let mut total_reference_minimizers: usize = 0usize;
        let mut total_seed_hits: usize = 0usize;
        let mut contig_names: Vec<ReferenceContigName> = Vec::new();
        files.reserve(reference_paths.len());
        contigs.reserve(reference_paths.len());

        if runtime_options.progress_enabled {
            emit_progress(
                "reference_build",
                &format!(
                    "event=start\tfiles_total={}\tsplit_n_run={split_n_run}",
                    reference_paths.len()
                ),
                build_start,
            );
        }
        check_memory_limit("reference build start", runtime_options)?;

        for (file_id, reference_path) in reference_paths.iter().enumerate() {
            let mut reader: fasta::io::Reader<Box<dyn io::BufRead>> =
                open_fasta_reader(reference_path)?;
            let mut mapped_length: u64 = 0u64;

            for result in reader.records() {
                let record: fasta::Record = result.map_err(|err| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!(
                            "failed to read FASTA record from reference {reference_path}: {err}"
                        ),
                    )
                })?;
                let record_name: String = String::from_utf8_lossy(record.name()).into_owned();
                let sequence: &fasta::record::Sequence = record.sequence();
                let sequence_bytes: &[u8] = sequence.as_ref();

                for segment_range in split_sequence_ranges(sequence_bytes, split_n_run) {
                    let segment_start: u32 = u32::try_from(segment_range.start).map_err(|err| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("reference segment start exceeds u32: {err}"),
                        )
                    })?;
                    let segment_end: u32 = u32::try_from(segment_range.end).map_err(|err| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("reference segment end exceeds u32: {err}"),
                        )
                    })?;
                    let segment_sequence: &[u8] = &sequence_bytes[segment_range];

                    mapped_length += mapped_length_from_fragment_ranges(
                        segment_sequence.len(),
                        fragment_length,
                        min_fragment_length,
                    );

                    let reference_contig_id: usize = contigs.len();
                    let sketch_capacity: usize =
                        segment_sequence.len() / (window_size / 2).max(1) + 1;
                    let mut reference_minimizers: Vec<ReferenceMinimizer> =
                        Vec::with_capacity(sketch_capacity);

                    if segment_sequence.len() >= kmer_size && segment_sequence.len() >= window_size
                    {
                        for (hash, position) in canonical_minimizers_with_positions(
                            segment_sequence,
                            kmer_size,
                            window_size,
                        ) {
                            reference_minimizers.push(ReferenceMinimizer { hash, position });
                        }
                    }

                    reference_minimizers.sort_unstable_by_key(|minimizer| minimizer.position);
                    index.reserve(reference_minimizers.len());

                    for minimizer in &reference_minimizers {
                        index.entry(minimizer.hash).or_default().push(SeedHit {
                            reference_contig_id: reference_contig_id as u32,
                            position: minimizer.position,
                        });
                    }
                    total_reference_minimizers += reference_minimizers.len();
                    total_seed_hits += reference_minimizers.len();

                    contigs.push(ReferenceContig {
                        file_id,
                        minimizers: reference_minimizers,
                    });
                    contig_names.push(ReferenceContigName {
                        file_id,
                        name: record_name.clone(),
                        segment_start,
                        segment_end,
                    });
                }
            }

            files.push(ReferenceFile {
                path: reference_path.clone(),
                mapped_length,
            });

            let files_done: usize = file_id + 1;
            if runtime_options.progress_enabled
                && (files_done % REFERENCE_PROGRESS_INTERVAL == 0
                    || files_done == reference_paths.len())
            {
                #[cfg(debug_assertions)]
                let estimated_struct_bytes: usize = reference_build_struct_bytes(
                    total_reference_minimizers,
                    total_seed_hits,
                    index.len(),
                );
                #[cfg(debug_assertions)]
                let progress_message: String = format!(
                    "event=files\tfiles_done={files_done}\tfiles_total={}\tcontigs={}\treference_minimizers={total_reference_minimizers}\tunique_minimizers={}\tseed_hits={total_seed_hits}\testimated_struct_mib={:.3}",
                    reference_paths.len(),
                    contigs.len(),
                    index.len(),
                    memory_mib(estimated_struct_bytes)
                );
                #[cfg(not(debug_assertions))]
                let progress_message: String = format!(
                    "event=files\tfiles_done={files_done}\tfiles_total={}\tcontigs={}\treference_minimizers={total_reference_minimizers}\tunique_minimizers={}\tseed_hits={total_seed_hits}",
                    reference_paths.len(),
                    contigs.len(),
                    index.len()
                );
                emit_progress("reference_build", &progress_message, build_start);
            }
            check_memory_limit(
                &format!(
                    "reference build after {files_done}/{} files",
                    reference_paths.len()
                ),
                runtime_options,
            )?;
        }

        if runtime_options.progress_enabled {
            emit_progress(
                "reference_build",
                &format!(
                    "event=complete\tfiles_done={}\tcontigs={}\treference_minimizers={total_reference_minimizers}\tunique_minimizers={}\tseed_hits={total_seed_hits}",
                    reference_paths.len(),
                    contigs.len(),
                    index.len()
                ),
                build_start,
            );
        }

        Ok(Self {
            files,
            contigs: ReferenceContigs::Owned(contigs),
            contig_names: Some(contig_names),
            index: ReferenceIndex::Hash(index),
        })
    }

    /// Build a reference sketch cache while streaming contig minimizers through scratch files.
    fn collect_and_save_streaming(
        reference_paths: &[String],
        kmer_size: usize,
        window_size: usize,
        fragment_length: u32,
        min_fragment_length: u32,
        split_n_run: usize,
        cache_path: &Path,
        tmp_dir: Option<&Path>,
        bgzip: bool,
        estimated_minimizers: usize,
        index_build_mode: IndexBuildMode,
        runtime_options: RuntimeOptions,
    ) -> io::Result<SketchBuildStats> {
        let effective_mode: IndexBuildMode =
            effective_index_build_mode(index_build_mode, estimated_minimizers);
        match effective_mode {
            IndexBuildMode::Hash => Self::collect_and_save_streaming_hash(
                reference_paths,
                kmer_size,
                window_size,
                fragment_length,
                min_fragment_length,
                split_n_run,
                cache_path,
                tmp_dir,
                bgzip,
                runtime_options,
            ),
            IndexBuildMode::Partitioned => Self::collect_and_save_streaming_partitioned(
                reference_paths,
                kmer_size,
                window_size,
                fragment_length,
                min_fragment_length,
                split_n_run,
                cache_path,
                tmp_dir,
                bgzip,
                estimated_minimizers,
                runtime_options,
            ),
            IndexBuildMode::Auto => unreachable!("auto mode is resolved before sketch build"),
        }
    }

    /// Build a reference sketch cache using the in-memory minimizer hit map.
    fn collect_and_save_streaming_hash(
        reference_paths: &[String],
        kmer_size: usize,
        window_size: usize,
        fragment_length: u32,
        min_fragment_length: u32,
        split_n_run: usize,
        cache_path: &Path,
        tmp_dir: Option<&Path>,
        bgzip: bool,
        runtime_options: RuntimeOptions,
    ) -> io::Result<SketchBuildStats> {
        let build_start: Instant = Instant::now();
        let mut files: Vec<ReferenceFile> = Vec::new();
        let mut contig_records: Vec<ContigRecord> = Vec::new();
        let mut contig_names: Vec<ReferenceContigName> = Vec::new();
        let mut index: ReferenceHitMap = ReferenceHitMap::default();
        let (reference_minimizer_scratch, reference_minimizer_file): (ScratchFile, fs::File) =
            ScratchFile::create(tmp_dir, "reference-build-minimizers")?;
        let mut reference_minimizer_writer: BufWriter<fs::File> =
            BufWriter::new(reference_minimizer_file);
        let mut reference_minimizer_count: usize = 0usize;
        let mut total_seed_hits: usize = 0usize;
        files.reserve(reference_paths.len());
        contig_records.reserve(reference_paths.len());

        if runtime_options.progress_enabled {
            emit_progress(
                "reference_build",
                &format!(
                    "event=start\tmode=streaming\tfiles_total={}\tsplit_n_run={split_n_run}\ttmp={}",
                    reference_paths.len(),
                    reference_minimizer_scratch.path.display()
                ),
                build_start,
            );
        }
        check_memory_limit("streaming reference build start", runtime_options)?;

        for (file_id, reference_path) in reference_paths.iter().enumerate() {
            let mut reader: fasta::io::Reader<Box<dyn io::BufRead>> =
                open_fasta_reader(reference_path)?;
            let mut mapped_length: u64 = 0u64;

            for result in reader.records() {
                let record: fasta::Record = result.map_err(|err| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!(
                            "failed to read FASTA record from reference {reference_path}: {err}"
                        ),
                    )
                })?;
                let record_name: String = String::from_utf8_lossy(record.name()).into_owned();
                let sequence: &fasta::record::Sequence = record.sequence();
                let sequence_bytes: &[u8] = sequence.as_ref();

                for segment_range in split_sequence_ranges(sequence_bytes, split_n_run) {
                    let segment_start: u32 = u32::try_from(segment_range.start).map_err(|err| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("reference segment start exceeds u32: {err}"),
                        )
                    })?;
                    let segment_end: u32 = u32::try_from(segment_range.end).map_err(|err| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("reference segment end exceeds u32: {err}"),
                        )
                    })?;
                    let segment_sequence: &[u8] = &sequence_bytes[segment_range];
                    mapped_length += mapped_length_from_fragment_ranges(
                        segment_sequence.len(),
                        fragment_length,
                        min_fragment_length,
                    );

                    let reference_contig_id: usize = contig_records.len();
                    let sketch_capacity: usize =
                        segment_sequence.len() / (window_size / 2).max(1) + 1;
                    let mut reference_minimizers: Vec<ReferenceMinimizer> =
                        Vec::with_capacity(sketch_capacity);

                    if segment_sequence.len() >= kmer_size && segment_sequence.len() >= window_size
                    {
                        for (hash, position) in canonical_minimizers_with_positions(
                            segment_sequence,
                            kmer_size,
                            window_size,
                        ) {
                            reference_minimizers.push(ReferenceMinimizer { hash, position });
                        }
                    }

                    reference_minimizers.sort_unstable_by_key(|minimizer| minimizer.position);
                    index.reserve(reference_minimizers.len());

                    for minimizer in &reference_minimizers {
                        index.entry(minimizer.hash).or_default().push(SeedHit {
                            reference_contig_id: reference_contig_id as u32,
                            position: minimizer.position,
                        });
                    }

                    let minimizer_offset: u64 = reference_minimizer_count as u64;
                    let minimizer_count: u32 =
                        u32::try_from(reference_minimizers.len()).map_err(|err| {
                            io::Error::new(
                                io::ErrorKind::InvalidData,
                                format!(
                                    "reference contig has too many minimizers for sketch cache: {err}"
                                ),
                            )
                        })?;
                    let file_id_u32: u32 = u32::try_from(file_id).map_err(|err| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("reference file id exceeds sketch cache limit: {err}"),
                        )
                    })?;

                    reference_minimizer_writer.write_all(slice_as_bytes(&reference_minimizers))?;
                    reference_minimizer_count += reference_minimizers.len();
                    total_seed_hits += reference_minimizers.len();
                    contig_records.push(ContigRecord {
                        minimizer_offset,
                        file_id: file_id_u32,
                        minimizer_count,
                    });
                    contig_names.push(ReferenceContigName {
                        file_id,
                        name: record_name.clone(),
                        segment_start,
                        segment_end,
                    });
                }
            }

            files.push(ReferenceFile {
                path: sketch_reference_name(reference_path),
                mapped_length,
            });

            let files_done: usize = file_id + 1;
            if runtime_options.progress_enabled
                && (files_done % REFERENCE_PROGRESS_INTERVAL == 0
                    || files_done == reference_paths.len())
            {
                emit_progress(
                    "reference_build",
                    &format!(
                        "event=files\tmode=streaming\tfiles_done={files_done}\tfiles_total={}\tcontigs={}\treference_minimizers={reference_minimizer_count}\tunique_minimizers={}\tseed_hits={total_seed_hits}",
                        reference_paths.len(),
                        contig_records.len(),
                        index.len()
                    ),
                    build_start,
                );
            }
            check_memory_limit(
                &format!(
                    "streaming reference build after {files_done}/{} files",
                    reference_paths.len()
                ),
                runtime_options,
            )?;
        }

        reference_minimizer_writer.flush()?;
        drop(reference_minimizer_writer);

        if runtime_options.progress_enabled {
            emit_progress(
                "reference_build",
                &format!(
                    "event=complete\tmode=streaming\tfiles_done={}\tcontigs={}\treference_minimizers={reference_minimizer_count}\tunique_minimizers={}\tseed_hits={total_seed_hits}",
                    reference_paths.len(),
                    contig_records.len(),
                    index.len()
                ),
                build_start,
            );
        }

        let stats: SketchBuildStats = SketchBuildStats {
            reference_count: files.len(),
            reference_contig_count: contig_records.len(),
            mapped_reference_length: files.iter().map(|file| file.mapped_length).sum(),
            reference_minimizer_count,
            unique_minimizer_count: index.len(),
        };

        Self::save_streamed_cache(
            cache_path,
            kmer_size,
            window_size,
            fragment_length,
            min_fragment_length,
            split_n_run,
            files,
            &index,
            contig_records,
            contig_names,
            reference_minimizer_count,
            &reference_minimizer_scratch,
            tmp_dir,
            bgzip,
            runtime_options,
        )?;

        Ok(stats)
    }

    fn save_streamed_cache(
        path: &Path,
        kmer_size: usize,
        window_size: usize,
        fragment_length: u32,
        min_fragment_length: u32,
        split_n_run: usize,
        files: Vec<ReferenceFile>,
        index: &ReferenceHitMap,
        contig_records: Vec<ContigRecord>,
        contig_names: Vec<ReferenceContigName>,
        reference_minimizer_count: usize,
        reference_minimizer_scratch: &ScratchFile,
        tmp_dir: Option<&Path>,
        bgzip: bool,
        runtime_options: RuntimeOptions,
    ) -> io::Result<()> {
        if let Some(parent) = path
            .parent()
            .filter(|parent| !parent.as_os_str().is_empty())
        {
            fs::create_dir_all(parent)?;
        }

        let save_start: Instant = Instant::now();
        if runtime_options.progress_enabled {
            emit_progress(
                "sketch_save",
                &format!(
                    "event=start\tmode=streaming\tunique_minimizers={}\tcontigs={}\tfiles={}\ttmp={}",
                    index.len(),
                    contig_records.len(),
                    files.len(),
                    reference_minimizer_scratch.path.display()
                ),
                save_start,
            );
        }
        check_memory_limit("streaming sketch save start", runtime_options)?;

        let keys: Vec<MinimizerKey> = index.keys().copied().collect::<Vec<_>>();
        let mphf: Mphf<MinimizerKey> = Mphf::new_parallel(1.7, &keys, None);
        let mut slot_keys: Vec<MinimizerKey> = vec![0; keys.len()];
        let mut hit_offsets: Vec<u64> = vec![0u64; keys.len()];
        let mut hit_counts: Vec<u32> = vec![0u32; keys.len()];
        let total_hits: usize = index.values().map(Vec::len).sum::<usize>();
        let mut hit_payloads: Vec<SeedHit> = Vec::with_capacity(total_hits);

        if runtime_options.progress_enabled {
            emit_progress(
                "sketch_save",
                &format!(
                    "event=arrays_allocated\tmode=streaming\tkey_count={}\thit_count={total_hits}\treference_minimizers={reference_minimizer_count}",
                    keys.len()
                ),
                save_start,
            );
        }
        check_memory_limit(
            "streaming sketch save after allocating arrays",
            runtime_options,
        )?;

        for (key_index, (key, hits)) in index.iter().enumerate() {
            let slot: usize = mphf.hash(key) as usize;
            slot_keys[slot] = *key;
            hit_offsets[slot] = hit_payloads.len() as u64;
            hit_counts[slot] = hits.len() as u32;
            hit_payloads.extend_from_slice(hits);

            let keys_done: usize = key_index + 1;
            if keys_done % SKETCH_KEY_PACK_PROGRESS_INTERVAL == 0 || keys_done == index.len() {
                if runtime_options.progress_enabled {
                    emit_progress(
                        "sketch_save",
                        &format!(
                            "event=pack_index\tmode=streaming\tkeys_done={keys_done}\tkey_count={}\thits_done={}",
                            index.len(),
                            hit_payloads.len()
                        ),
                        save_start,
                    );
                }
                check_memory_limit("streaming sketch save while packing index", runtime_options)?;
            }
        }

        let expected_reference_minimizer_bytes: u64 =
            (reference_minimizer_count * size_of::<ReferenceMinimizer>()) as u64;
        let actual_reference_minimizer_bytes: u64 =
            fs::metadata(&reference_minimizer_scratch.path)?.len();
        if actual_reference_minimizer_bytes != expected_reference_minimizer_bytes {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "reference minimizer scratch file has {actual_reference_minimizer_bytes} bytes, expected {expected_reference_minimizer_bytes}"
                ),
            ));
        }
        if contig_names.len() != contig_records.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "contig sidecar would have {} records, expected {}",
                    contig_names.len(),
                    contig_records.len()
                ),
            ));
        }
        write_contig_name_sidecar(
            path,
            &files,
            &contig_names,
            runtime_options.effective_worker_threads(),
        )?;

        let metadata: CachedReferenceMetadata = CachedReferenceMetadata {
            version: SKETCH_VERSION,
            k: kmer_size,
            w: window_size,
            key_mode: SKETCH_KEY_MODE.to_string(),
            fragment_length,
            min_fragment_length,
            split_n_run,
            dust_enabled: false,
            files,
            mphf,
            key_count: slot_keys.len(),
            hit_count: hit_payloads.len(),
            contig_count: contig_records.len(),
            reference_minimizer_count,
        };
        let metadata_bytes: Vec<u8> = serde_json::to_vec(&metadata).map_err(|err| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("failed to encode reference sketch metadata: {err}"),
            )
        })?;

        let metadata_end: usize = SKETCH_MAGIC
            .len()
            .checked_add(size_of::<u64>())
            .and_then(|offset| offset.checked_add(metadata_bytes.len()))
            .ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "metadata length overflow")
            })?;
        let slot_keys_offset: usize = align_up(metadata_end, 8);
        let hit_offsets_offset: usize = align_up(
            checked_section_end(slot_keys_offset, slot_keys.len(), size_of::<MinimizerKey>())?,
            align_of::<u64>(),
        );
        let hit_counts_offset: usize =
            checked_section_end(hit_offsets_offset, hit_offsets.len(), size_of::<u64>())?;
        let hit_payloads_offset: usize = align_up(
            checked_section_end(hit_counts_offset, hit_counts.len(), size_of::<u32>())?,
            align_of::<SeedHit>(),
        );
        let contig_records_offset: usize = align_up(
            checked_section_end(
                hit_payloads_offset,
                hit_payloads.len(),
                size_of::<SeedHit>(),
            )?,
            align_of::<ContigRecord>(),
        );
        let reference_minimizers_offset: usize = align_up(
            checked_section_end(
                contig_records_offset,
                contig_records.len(),
                size_of::<ContigRecord>(),
            )?,
            align_of::<ReferenceMinimizer>(),
        );

        let mut sketch_output: SketchOutput = SketchOutput::create(
            path,
            tmp_dir,
            bgzip,
            runtime_options.effective_worker_threads(),
        )?;
        let mut writer: &mut BufWriter<fs::File> = sketch_output.writer_mut()?;
        writer.write_all(SKETCH_MAGIC)?;
        writer.write_all(&(metadata_bytes.len() as u64).to_le_bytes())?;
        writer.write_all(&metadata_bytes)?;
        write_padding(&mut writer, slot_keys_offset - metadata_end)?;
        writer.write_all(slice_as_bytes(&slot_keys))?;
        write_padding(
            &mut writer,
            hit_offsets_offset
                - checked_section_end(
                    slot_keys_offset,
                    slot_keys.len(),
                    size_of::<MinimizerKey>(),
                )?,
        )?;
        writer.write_all(slice_as_bytes(&hit_offsets))?;
        writer.write_all(slice_as_bytes(&hit_counts))?;
        write_padding(
            &mut writer,
            hit_payloads_offset
                - checked_section_end(hit_counts_offset, hit_counts.len(), size_of::<u32>())?,
        )?;
        writer.write_all(slice_as_bytes(&hit_payloads))?;
        write_padding(
            &mut writer,
            contig_records_offset
                - checked_section_end(
                    hit_payloads_offset,
                    hit_payloads.len(),
                    size_of::<SeedHit>(),
                )?,
        )?;
        writer.write_all(slice_as_bytes(&contig_records))?;
        write_padding(
            &mut writer,
            reference_minimizers_offset
                - checked_section_end(
                    contig_records_offset,
                    contig_records.len(),
                    size_of::<ContigRecord>(),
                )?,
        )?;
        let mut reference_minimizer_reader: fs::File =
            fs::File::open(&reference_minimizer_scratch.path)?;
        io::copy(&mut reference_minimizer_reader, &mut writer)?;
        let output_file_bytes: u64 = sketch_output.finish()?;

        if runtime_options.progress_enabled {
            emit_progress(
                "sketch_save",
                &format!(
                    "event=complete\tmode=streaming\tpath={}\tfile_bytes={}",
                    path.display(),
                    output_file_bytes
                ),
                save_start,
            );
        }
        check_memory_limit("streaming sketch save complete", runtime_options)
    }

    /// Build a reference sketch cache using disk-partitioned minimizer hit records.
    fn collect_and_save_streaming_partitioned(
        reference_paths: &[String],
        kmer_size: usize,
        window_size: usize,
        fragment_length: u32,
        min_fragment_length: u32,
        split_n_run: usize,
        cache_path: &Path,
        tmp_dir: Option<&Path>,
        bgzip: bool,
        estimated_minimizers: usize,
        runtime_options: RuntimeOptions,
    ) -> io::Result<SketchBuildStats> {
        let build_start: Instant = Instant::now();
        let partition_plan: PartitionBuildPlan =
            partition_build_plan(estimated_minimizers, runtime_options.max_memory_bytes);
        let mut files: Vec<ReferenceFile> = Vec::new();
        let mut contig_records: Vec<ContigRecord> = Vec::new();
        let mut contig_names: Vec<ReferenceContigName> = Vec::new();
        let mut partition_writers: PartitionWriters = PartitionWriters::new(
            partition_plan.partition_count,
            tmp_dir,
            PARTITION_BUFFER_RECORDS,
        )?;
        let (reference_minimizer_scratch, reference_minimizer_file): (ScratchFile, fs::File) =
            ScratchFile::create(tmp_dir, "reference-build-minimizers")?;
        let mut reference_minimizer_writer: BufWriter<fs::File> =
            BufWriter::new(reference_minimizer_file);
        let mut reference_minimizer_count: usize = 0usize;
        let mut total_seed_hits: usize = 0usize;
        files.reserve(reference_paths.len());
        contig_records.reserve(reference_paths.len());

        if runtime_options.progress_enabled {
            emit_progress(
                "reference_build",
                &format!(
                    "event=start\tmode=streaming\tindex_build_mode=partitioned\tfiles_total={}\tsplit_n_run={split_n_run}\tpartitions={}\testimated_minimizers={estimated_minimizers}\testimated_record_mib={:.3}\ttarget_partition_mib={:.3}\ttmp={}",
                    reference_paths.len(),
                    partition_plan.partition_count,
                    memory_mib(partition_plan.estimated_record_bytes),
                    memory_mib(partition_plan.target_partition_bytes),
                    reference_minimizer_scratch.path.display()
                ),
                build_start,
            );
        }
        check_memory_limit("partitioned reference build start", runtime_options)?;

        for (file_id, reference_path) in reference_paths.iter().enumerate() {
            let mut reader: fasta::io::Reader<Box<dyn io::BufRead>> =
                open_fasta_reader(reference_path)?;
            let mut mapped_length: u64 = 0u64;

            for result in reader.records() {
                let record: fasta::Record = result.map_err(|err| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!(
                            "failed to read FASTA record from reference {reference_path}: {err}"
                        ),
                    )
                })?;
                let record_name: String = String::from_utf8_lossy(record.name()).into_owned();
                let sequence: &fasta::record::Sequence = record.sequence();
                let sequence_bytes: &[u8] = sequence.as_ref();

                for segment_range in split_sequence_ranges(sequence_bytes, split_n_run) {
                    let segment_start: u32 = u32::try_from(segment_range.start).map_err(|err| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("reference segment start exceeds u32: {err}"),
                        )
                    })?;
                    let segment_end: u32 = u32::try_from(segment_range.end).map_err(|err| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("reference segment end exceeds u32: {err}"),
                        )
                    })?;
                    let segment_sequence: &[u8] = &sequence_bytes[segment_range];
                    mapped_length += mapped_length_from_fragment_ranges(
                        segment_sequence.len(),
                        fragment_length,
                        min_fragment_length,
                    );

                    let reference_contig_id: usize = contig_records.len();
                    let sketch_capacity: usize =
                        segment_sequence.len() / (window_size / 2).max(1) + 1;
                    let mut reference_minimizers: Vec<ReferenceMinimizer> =
                        Vec::with_capacity(sketch_capacity);

                    if segment_sequence.len() >= kmer_size && segment_sequence.len() >= window_size
                    {
                        for (hash, position) in canonical_minimizers_with_positions(
                            segment_sequence,
                            kmer_size,
                            window_size,
                        ) {
                            reference_minimizers.push(ReferenceMinimizer { hash, position });
                        }
                    }

                    reference_minimizers.sort_unstable_by_key(|minimizer| minimizer.position);
                    let reference_contig_id_u32: u32 =
                        u32::try_from(reference_contig_id).map_err(|err| {
                            io::Error::new(
                                io::ErrorKind::InvalidData,
                                format!("reference contig id exceeds sketch cache limit: {err}"),
                            )
                        })?;
                    for minimizer in &reference_minimizers {
                        partition_writers.push(PartitionHitRecord {
                            key: minimizer.hash,
                            hit: SeedHit {
                                reference_contig_id: reference_contig_id_u32,
                                position: minimizer.position,
                            },
                        })?;
                    }

                    let minimizer_offset: u64 = reference_minimizer_count as u64;
                    let minimizer_count: u32 =
                        u32::try_from(reference_minimizers.len()).map_err(|err| {
                            io::Error::new(
                                io::ErrorKind::InvalidData,
                                format!(
                                    "reference contig has too many minimizers for sketch cache: {err}"
                                ),
                            )
                        })?;
                    let file_id_u32: u32 = u32::try_from(file_id).map_err(|err| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("reference file id exceeds sketch cache limit: {err}"),
                        )
                    })?;

                    reference_minimizer_writer.write_all(slice_as_bytes(&reference_minimizers))?;
                    reference_minimizer_count += reference_minimizers.len();
                    total_seed_hits += reference_minimizers.len();
                    contig_records.push(ContigRecord {
                        minimizer_offset,
                        file_id: file_id_u32,
                        minimizer_count,
                    });
                    contig_names.push(ReferenceContigName {
                        file_id,
                        name: record_name.clone(),
                        segment_start,
                        segment_end,
                    });
                }
            }

            files.push(ReferenceFile {
                path: sketch_reference_name(reference_path),
                mapped_length,
            });

            let files_done: usize = file_id + 1;
            if runtime_options.progress_enabled
                && (files_done % REFERENCE_PROGRESS_INTERVAL == 0
                    || files_done == reference_paths.len())
            {
                emit_progress(
                    "reference_build",
                    &format!(
                        "event=files\tmode=streaming\tindex_build_mode=partitioned\tfiles_done={files_done}\tfiles_total={}\tcontigs={}\treference_minimizers={reference_minimizer_count}\tseed_hits={total_seed_hits}\tpartitions={}",
                        reference_paths.len(),
                        contig_records.len(),
                        partition_plan.partition_count
                    ),
                    build_start,
                );
            }
            check_memory_limit(
                &format!(
                    "partitioned reference build after {files_done}/{} files",
                    reference_paths.len()
                ),
                runtime_options,
            )?;
        }

        reference_minimizer_writer.flush()?;
        drop(reference_minimizer_writer);
        partition_writers.flush_all()?;

        if runtime_options.progress_enabled {
            emit_progress(
                "reference_build",
                &format!(
                    "event=complete\tmode=streaming\tindex_build_mode=partitioned\tfiles_done={}\tcontigs={}\treference_minimizers={reference_minimizer_count}\tseed_hits={total_seed_hits}\tpartitions={}",
                    reference_paths.len(),
                    contig_records.len(),
                    partition_plan.partition_count
                ),
                build_start,
            );
        }

        let reference_count: usize = files.len();
        let reference_contig_count: usize = contig_records.len();
        let mapped_reference_length: u64 = files.iter().map(|file| file.mapped_length).sum();
        let unique_minimizer_count: usize = Self::save_partitioned_streamed_cache(
            cache_path,
            kmer_size,
            window_size,
            fragment_length,
            min_fragment_length,
            split_n_run,
            files,
            contig_records,
            contig_names,
            reference_minimizer_count,
            &reference_minimizer_scratch,
            &partition_writers,
            partition_plan,
            tmp_dir,
            bgzip,
            runtime_options,
        )?;

        Ok(SketchBuildStats {
            reference_count,
            reference_contig_count,
            mapped_reference_length,
            reference_minimizer_count,
            unique_minimizer_count,
        })
    }

    fn read_partition_records(path: &Path) -> io::Result<Vec<PartitionHitRecord>> {
        let byte_len: u64 = fs::metadata(path)?.len();
        let record_size: u64 = size_of::<PartitionHitRecord>() as u64;
        if byte_len % record_size != 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "partition file {} has {byte_len} bytes, not a multiple of {record_size}",
                    path.display()
                ),
            ));
        }

        let record_count: usize = usize::try_from(byte_len / record_size).map_err(|err| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("partition file has too many records to load: {err}"),
            )
        })?;
        let mut records: Vec<PartitionHitRecord> =
            vec![PartitionHitRecord::default(); record_count];
        if record_count > 0 {
            let mut reader: BufReader<fs::File> = BufReader::new(fs::File::open(path)?);
            reader.read_exact(slice_as_bytes_mut(&mut records))?;
        }

        Ok(records)
    }

    fn group_partition_records(
        partition_index: usize,
        partition_path: &Path,
        tmp_dir: Option<&Path>,
    ) -> io::Result<PartitionGroupResult> {
        let mut records: Vec<PartitionHitRecord> = Self::read_partition_records(partition_path)?;
        records.sort_unstable_by_key(|record| {
            (
                record.key,
                record.hit.reference_contig_id,
                record.hit.position,
            )
        });

        let (grouped_key_scratch, grouped_key_file): (ScratchFile, fs::File) = ScratchFile::create(
            tmp_dir,
            &format!("partition-{partition_index}-grouped-keys"),
        )?;
        let (hit_payload_scratch, hit_payload_file): (ScratchFile, fs::File) = ScratchFile::create(
            tmp_dir,
            &format!("partition-{partition_index}-hit-payloads"),
        )?;
        let mut grouped_key_writer: BufWriter<fs::File> = BufWriter::new(grouped_key_file);
        let mut hit_payload_writer: BufWriter<fs::File> = BufWriter::new(hit_payload_file);
        let mut keys: Vec<MinimizerKey> = Vec::new();
        let mut local_hit_count: usize = 0usize;
        let mut hit_buffer: Vec<SeedHit> = Vec::with_capacity(1_048_576);

        let mut group_start: usize = 0usize;
        while group_start < records.len() {
            let key: MinimizerKey = records[group_start].key;
            let mut group_end: usize = group_start + 1;
            while group_end < records.len() && records[group_end].key == key {
                group_end += 1;
            }

            let hit_count: usize = group_end - group_start;
            let hit_count_u32: u32 = u32::try_from(hit_count).map_err(|err| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("minimizer hit count exceeds sketch cache limit: {err}"),
                )
            })?;
            let hit_offset_u64: u64 = u64::try_from(local_hit_count).map_err(|err| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("minimizer hit offset exceeds sketch cache limit: {err}"),
                )
            })?;
            let grouped_record: GroupedKeyRecord = GroupedKeyRecord {
                key,
                hit_offset: hit_offset_u64,
                hit_count: hit_count_u32,
            };
            grouped_key_writer.write_all(slice_as_bytes(std::slice::from_ref(&grouped_record)))?;
            keys.push(key);

            for record in &records[group_start..group_end] {
                hit_buffer.push(record.hit);
                if hit_buffer.len() >= 1_048_576 {
                    hit_payload_writer.write_all(slice_as_bytes(&hit_buffer))?;
                    hit_buffer.clear();
                }
            }

            local_hit_count = local_hit_count.checked_add(hit_count).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "partition hit count overflow")
            })?;
            group_start = group_end;
        }

        if !hit_buffer.is_empty() {
            hit_payload_writer.write_all(slice_as_bytes(&hit_buffer))?;
        }
        grouped_key_writer.flush()?;
        hit_payload_writer.flush()?;
        drop(grouped_key_writer);
        drop(hit_payload_writer);

        let expected_grouped_key_bytes: u64 = (keys.len() * size_of::<GroupedKeyRecord>()) as u64;
        let actual_grouped_key_bytes: u64 = fs::metadata(&grouped_key_scratch.path)?.len();
        if actual_grouped_key_bytes != expected_grouped_key_bytes {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "partition grouped key scratch file has {actual_grouped_key_bytes} bytes, expected {expected_grouped_key_bytes}"
                ),
            ));
        }

        let expected_hit_payload_bytes: u64 = (local_hit_count * size_of::<SeedHit>()) as u64;
        let actual_hit_payload_bytes: u64 = fs::metadata(&hit_payload_scratch.path)?.len();
        if actual_hit_payload_bytes != expected_hit_payload_bytes {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "partition hit payload scratch file has {actual_hit_payload_bytes} bytes, expected {expected_hit_payload_bytes}"
                ),
            ));
        }

        Ok(PartitionGroupResult {
            partition_index,
            grouped_key_scratch,
            hit_payload_scratch,
            keys,
            hit_count: local_hit_count,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn save_partitioned_streamed_cache(
        path: &Path,
        kmer_size: usize,
        window_size: usize,
        fragment_length: u32,
        min_fragment_length: u32,
        split_n_run: usize,
        files: Vec<ReferenceFile>,
        contig_records: Vec<ContigRecord>,
        contig_names: Vec<ReferenceContigName>,
        reference_minimizer_count: usize,
        reference_minimizer_scratch: &ScratchFile,
        partition_writers: &PartitionWriters,
        partition_plan: PartitionBuildPlan,
        tmp_dir: Option<&Path>,
        bgzip: bool,
        runtime_options: RuntimeOptions,
    ) -> io::Result<usize> {
        if let Some(parent) = path
            .parent()
            .filter(|parent| !parent.as_os_str().is_empty())
        {
            fs::create_dir_all(parent)?;
        }

        let save_start: Instant = Instant::now();
        if runtime_options.progress_enabled {
            emit_progress(
                "sketch_save",
                &format!(
                    "event=start\tmode=streaming\tindex_build_mode=partitioned\tpartitions={}\testimated_minimizers={}\ttarget_partition_mib={:.3}\tcontigs={}\tfiles={}\ttmp={}",
                    partition_plan.partition_count,
                    partition_plan.estimated_record_bytes / size_of::<PartitionHitRecord>(),
                    memory_mib(partition_plan.target_partition_bytes),
                    contig_records.len(),
                    files.len(),
                    reference_minimizer_scratch.path.display()
                ),
                save_start,
            );
        }
        check_memory_limit("partitioned sketch save start", runtime_options)?;

        let partition_paths: Vec<PathBuf> = (0..partition_writers.partition_count())
            .map(|partition_index| partition_writers.path(partition_index).to_path_buf())
            .collect();
        let partition_parallelism: usize = runtime_options
            .effective_worker_threads()
            .min(partition_paths.len().max(1));
        let completed_partitions: AtomicUsize = AtomicUsize::new(0);
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(partition_parallelism)
            .build()
            .map_err(|err| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("failed to initialize partition grouping thread pool: {err}"),
                )
            })?;
        let mut partition_results: Vec<PartitionGroupResult> = pool.install(|| {
            partition_paths
                .par_iter()
                .enumerate()
                .map(|(partition_index, partition_path)| {
                    let result: PartitionGroupResult =
                        Self::group_partition_records(partition_index, partition_path, tmp_dir)?;
                    let partitions_done: usize =
                        completed_partitions.fetch_add(1, AtomicOrdering::Relaxed) + 1;
                    if runtime_options.progress_enabled
                        && (partitions_done % 16 == 0 || partitions_done == partition_paths.len())
                    {
                        emit_progress(
                            "sketch_save",
                            &format!(
                                "event=sort_group_partitions\tmode=streaming\tindex_build_mode=partitioned\tpartitions_done={partitions_done}\tpartitions={}\tpartition_parallelism={partition_parallelism}",
                                partition_paths.len()
                            ),
                            save_start,
                        );
                    }
                    check_memory_limit(
                        "partitioned sketch save while sorting partitions",
                        runtime_options,
                    )?;

                    Ok(result)
                })
                .collect::<io::Result<Vec<_>>>()
        })?;
        partition_results.sort_by_key(|result| result.partition_index);

        let total_unique_minimizers: usize = partition_results
            .iter()
            .map(|result| result.keys.len())
            .sum();
        let mut keys: Vec<MinimizerKey> = Vec::with_capacity(total_unique_minimizers);
        let mut total_hits: usize = 0usize;
        for result in &partition_results {
            keys.extend_from_slice(&result.keys);
            total_hits = total_hits.checked_add(result.hit_count).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "total hit count overflow")
            })?;
        }

        if runtime_options.progress_enabled {
            emit_progress(
                "sketch_save",
                &format!(
                    "event=keys_collected\tmode=streaming\tindex_build_mode=partitioned\tkey_count={}\thit_count={total_hits}\tpartition_parallelism={partition_parallelism}",
                    keys.len(),
                ),
                save_start,
            );
        }
        check_memory_limit(
            "partitioned sketch save after collecting keys",
            runtime_options,
        )?;

        let mphf: Mphf<MinimizerKey> = Mphf::new_parallel(1.7, &keys, None);
        if runtime_options.progress_enabled {
            emit_progress(
                "sketch_save",
                &format!(
                    "event=mphf_built\tmode=streaming\tindex_build_mode=partitioned\tkey_count={}",
                    keys.len()
                ),
                save_start,
            );
        }
        check_memory_limit(
            "partitioned sketch save after building MPH",
            runtime_options,
        )?;

        let mut slot_keys: Vec<MinimizerKey> = vec![0; keys.len()];
        let mut hit_offsets: Vec<u64> = vec![0u64; keys.len()];
        let mut hit_counts: Vec<u32> = vec![0u32; keys.len()];
        let mut grouped_records_done: usize = 0usize;
        let mut partition_hit_offset: u64 = 0u64;
        for result in &partition_results {
            let mut grouped_key_reader: BufReader<fs::File> =
                BufReader::new(fs::File::open(&result.grouped_key_scratch.path)?);
            let mut partition_records_done: usize = 0usize;
            while partition_records_done < result.keys.len() {
                let records_to_read: usize =
                    (result.keys.len() - partition_records_done).min(1_000_000);
                let mut grouped_records: Vec<GroupedKeyRecord> =
                    vec![GroupedKeyRecord::default(); records_to_read];
                grouped_key_reader.read_exact(slice_as_bytes_mut(&mut grouped_records))?;

                for grouped_record in grouped_records {
                    let slot: usize = mphf.hash(&grouped_record.key) as usize;
                    slot_keys[slot] = grouped_record.key;
                    hit_offsets[slot] = partition_hit_offset
                        .checked_add(grouped_record.hit_offset)
                        .ok_or_else(|| {
                            io::Error::new(io::ErrorKind::InvalidData, "global hit offset overflow")
                        })?;
                    hit_counts[slot] = grouped_record.hit_count;
                }

                partition_records_done += records_to_read;
                grouped_records_done += records_to_read;
                if runtime_options.progress_enabled
                    && (grouped_records_done % SKETCH_KEY_PACK_PROGRESS_INTERVAL == 0
                        || grouped_records_done == keys.len())
                {
                    emit_progress(
                        "sketch_save",
                        &format!(
                            "event=pack_index\tmode=streaming\tindex_build_mode=partitioned\tkeys_done={grouped_records_done}\tkey_count={}",
                            keys.len()
                        ),
                        save_start,
                    );
                }
                check_memory_limit(
                    "partitioned sketch save while packing index",
                    runtime_options,
                )?;
            }
            partition_hit_offset = partition_hit_offset
                .checked_add(u64::try_from(result.hit_count).map_err(|err| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("partition hit count exceeds u64: {err}"),
                    )
                })?)
                .ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidData, "partition hit offset overflow")
                })?;
        }

        let expected_reference_minimizer_bytes: u64 =
            (reference_minimizer_count * size_of::<ReferenceMinimizer>()) as u64;
        let actual_reference_minimizer_bytes: u64 =
            fs::metadata(&reference_minimizer_scratch.path)?.len();
        if actual_reference_minimizer_bytes != expected_reference_minimizer_bytes {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "reference minimizer scratch file has {actual_reference_minimizer_bytes} bytes, expected {expected_reference_minimizer_bytes}"
                ),
            ));
        }
        if contig_names.len() != contig_records.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "contig sidecar would have {} records, expected {}",
                    contig_names.len(),
                    contig_records.len()
                ),
            ));
        }
        write_contig_name_sidecar(
            path,
            &files,
            &contig_names,
            runtime_options.effective_worker_threads(),
        )?;

        let metadata: CachedReferenceMetadata = CachedReferenceMetadata {
            version: SKETCH_VERSION,
            k: kmer_size,
            w: window_size,
            key_mode: SKETCH_KEY_MODE.to_string(),
            fragment_length,
            min_fragment_length,
            split_n_run,
            dust_enabled: false,
            files,
            mphf,
            key_count: slot_keys.len(),
            hit_count: total_hits,
            contig_count: contig_records.len(),
            reference_minimizer_count,
        };
        let metadata_bytes: Vec<u8> = serde_json::to_vec(&metadata).map_err(|err| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("failed to encode reference sketch metadata: {err}"),
            )
        })?;

        let metadata_end: usize = SKETCH_MAGIC
            .len()
            .checked_add(size_of::<u64>())
            .and_then(|offset| offset.checked_add(metadata_bytes.len()))
            .ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "metadata length overflow")
            })?;
        let slot_keys_offset: usize = align_up(metadata_end, 8);
        let hit_offsets_offset: usize = align_up(
            checked_section_end(slot_keys_offset, slot_keys.len(), size_of::<MinimizerKey>())?,
            align_of::<u64>(),
        );
        let hit_counts_offset: usize =
            checked_section_end(hit_offsets_offset, hit_offsets.len(), size_of::<u64>())?;
        let hit_payloads_offset: usize = align_up(
            checked_section_end(hit_counts_offset, hit_counts.len(), size_of::<u32>())?,
            align_of::<SeedHit>(),
        );
        let contig_records_offset: usize = align_up(
            checked_section_end(hit_payloads_offset, total_hits, size_of::<SeedHit>())?,
            align_of::<ContigRecord>(),
        );
        let reference_minimizers_offset: usize = align_up(
            checked_section_end(
                contig_records_offset,
                contig_records.len(),
                size_of::<ContigRecord>(),
            )?,
            align_of::<ReferenceMinimizer>(),
        );

        let mut sketch_output: SketchOutput = SketchOutput::create(
            path,
            tmp_dir,
            bgzip,
            runtime_options.effective_worker_threads(),
        )?;
        let mut writer: &mut BufWriter<fs::File> = sketch_output.writer_mut()?;
        writer.write_all(SKETCH_MAGIC)?;
        writer.write_all(&(metadata_bytes.len() as u64).to_le_bytes())?;
        writer.write_all(&metadata_bytes)?;
        write_padding(&mut writer, slot_keys_offset - metadata_end)?;
        writer.write_all(slice_as_bytes(&slot_keys))?;
        write_padding(
            &mut writer,
            hit_offsets_offset
                - checked_section_end(
                    slot_keys_offset,
                    slot_keys.len(),
                    size_of::<MinimizerKey>(),
                )?,
        )?;
        writer.write_all(slice_as_bytes(&hit_offsets))?;
        writer.write_all(slice_as_bytes(&hit_counts))?;
        write_padding(
            &mut writer,
            hit_payloads_offset
                - checked_section_end(hit_counts_offset, hit_counts.len(), size_of::<u32>())?,
        )?;
        for result in &partition_results {
            let mut hit_payload_reader: fs::File =
                fs::File::open(&result.hit_payload_scratch.path)?;
            io::copy(&mut hit_payload_reader, &mut writer)?;
        }
        write_padding(
            &mut writer,
            contig_records_offset
                - checked_section_end(hit_payloads_offset, total_hits, size_of::<SeedHit>())?,
        )?;
        writer.write_all(slice_as_bytes(&contig_records))?;
        write_padding(
            &mut writer,
            reference_minimizers_offset
                - checked_section_end(
                    contig_records_offset,
                    contig_records.len(),
                    size_of::<ContigRecord>(),
                )?,
        )?;
        let mut reference_minimizer_reader: fs::File =
            fs::File::open(&reference_minimizer_scratch.path)?;
        io::copy(&mut reference_minimizer_reader, &mut writer)?;
        let output_file_bytes: u64 = sketch_output.finish()?;

        if runtime_options.progress_enabled {
            emit_progress(
                "sketch_save",
                &format!(
                    "event=complete\tmode=streaming\tindex_build_mode=partitioned\tpath={}\tfile_bytes={}",
                    path.display(),
                    output_file_bytes
                ),
                save_start,
            );
        }
        check_memory_limit("partitioned sketch save complete", runtime_options)?;

        Ok(keys.len())
    }

    /// Save the in-memory reference index as a zero-copy-loadable sketch cache.
    #[cfg(test)]
    fn save(
        &self,
        path: &Path,
        kmer_size: usize,
        window_size: usize,
        fragment_length: u32,
        min_fragment_length: u32,
        split_n_run: usize,
        tmp_dir: Option<&Path>,
        runtime_options: RuntimeOptions,
    ) -> io::Result<()> {
        let ReferenceIndex::Hash(index) = &self.index else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "reference sketch cache can only be saved from an in-memory HashMap index",
            ));
        };
        let ReferenceContigs::Owned(contigs) = &self.contigs else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "reference sketch cache can only be saved from owned reference contigs",
            ));
        };

        if let Some(parent) = path
            .parent()
            .filter(|parent| !parent.as_os_str().is_empty())
        {
            fs::create_dir_all(parent)?;
        }

        let save_start: Instant = Instant::now();
        if runtime_options.progress_enabled {
            emit_progress(
                "sketch_save",
                &format!(
                    "event=start\tunique_minimizers={}\tcontigs={}\tfiles={}",
                    index.len(),
                    contigs.len(),
                    self.files.len()
                ),
                save_start,
            );
        }
        check_memory_limit("sketch save start", runtime_options)?;

        let keys: Vec<MinimizerKey> = index.keys().copied().collect::<Vec<_>>();
        if runtime_options.progress_enabled {
            emit_progress(
                "sketch_save",
                &format!("event=keys_collected\tkey_count={}", keys.len()),
                save_start,
            );
        }
        check_memory_limit("sketch save after collecting keys", runtime_options)?;
        let mphf: Mphf<MinimizerKey> = Mphf::new_parallel(1.7, &keys, None);
        if runtime_options.progress_enabled {
            emit_progress(
                "sketch_save",
                &format!("event=mphf_built\tkey_count={}", keys.len()),
                save_start,
            );
        }
        check_memory_limit("sketch save after building MPH", runtime_options)?;
        let mut slot_keys: Vec<MinimizerKey> = vec![0; keys.len()];
        let mut hit_offsets: Vec<u64> = vec![0u64; keys.len()];
        let mut hit_counts: Vec<u32> = vec![0u32; keys.len()];
        let total_hits: usize = index.values().map(Vec::len).sum::<usize>();
        let mut hit_payloads: Vec<SeedHit> = Vec::with_capacity(total_hits);
        let total_reference_minimizers: usize =
            contigs.iter().map(|contig| contig.minimizers.len()).sum();
        let mut contig_records: Vec<ContigRecord> = Vec::with_capacity(contigs.len());
        let (reference_minimizer_scratch, reference_minimizer_file): (ScratchFile, fs::File) =
            ScratchFile::create(tmp_dir, "reference-minimizers")?;
        let mut reference_minimizer_writer: BufWriter<fs::File> =
            BufWriter::new(reference_minimizer_file);
        let mut reference_minimizer_count: usize = 0usize;
        let files: Vec<ReferenceFile> = self
            .files
            .iter()
            .map(|file| ReferenceFile {
                path: sketch_reference_name(&file.path),
                mapped_length: file.mapped_length,
            })
            .collect();

        if runtime_options.progress_enabled {
            #[cfg(debug_assertions)]
            let estimated_pack_bytes: usize = keys
                .len()
                .saturating_mul(size_of::<MinimizerKey>() + size_of::<u64>() + size_of::<u32>())
                .saturating_add(total_hits.saturating_mul(size_of::<SeedHit>()))
                .saturating_add(contigs.len().saturating_mul(size_of::<ContigRecord>()))
                .saturating_add(
                    total_reference_minimizers.saturating_mul(size_of::<ReferenceMinimizer>()),
                );
            #[cfg(debug_assertions)]
            let progress_message: String = format!(
                "event=arrays_allocated\tkey_count={}\thit_count={total_hits}\treference_minimizers={total_reference_minimizers}\testimated_pack_mib={:.3}\ttmp={}",
                keys.len(),
                memory_mib(estimated_pack_bytes),
                reference_minimizer_scratch.path.display()
            );
            #[cfg(not(debug_assertions))]
            let progress_message: String = format!(
                "event=arrays_allocated\tkey_count={}\thit_count={total_hits}\treference_minimizers={total_reference_minimizers}\ttmp={}",
                keys.len(),
                reference_minimizer_scratch.path.display()
            );
            emit_progress("sketch_save", &progress_message, save_start);
        }
        check_memory_limit("sketch save after allocating arrays", runtime_options)?;

        for (key_index, (key, hits)) in index.iter().enumerate() {
            let slot: usize = mphf.hash(key) as usize;
            slot_keys[slot] = *key;
            hit_offsets[slot] = hit_payloads.len() as u64;
            hit_counts[slot] = hits.len() as u32;
            hit_payloads.extend_from_slice(hits);

            let keys_done: usize = key_index + 1;
            if keys_done % SKETCH_KEY_PACK_PROGRESS_INTERVAL == 0 || keys_done == index.len() {
                if runtime_options.progress_enabled {
                    emit_progress(
                        "sketch_save",
                        &format!(
                            "event=pack_index\tkeys_done={keys_done}\tkey_count={}\thits_done={}",
                            index.len(),
                            hit_payloads.len()
                        ),
                        save_start,
                    );
                }
                check_memory_limit("sketch save while packing index", runtime_options)?;
            }
        }

        for (contig_index, contig) in contigs.iter().enumerate() {
            let minimizer_offset: u64 = reference_minimizer_count as u64;
            let minimizer_count: u32 = u32::try_from(contig.minimizers.len()).map_err(|err| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("reference contig has too many minimizers for sketch cache: {err}"),
                )
            })?;
            let file_id: u32 = u32::try_from(contig.file_id).map_err(|err| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("reference file id exceeds sketch cache limit: {err}"),
                )
            })?;

            contig_records.push(ContigRecord {
                minimizer_offset,
                file_id,
                minimizer_count,
            });
            reference_minimizer_writer.write_all(slice_as_bytes(&contig.minimizers))?;
            reference_minimizer_count += contig.minimizers.len();

            let contigs_done: usize = contig_index + 1;
            if contigs_done % SKETCH_CONTIG_PACK_PROGRESS_INTERVAL == 0
                || contigs_done == contigs.len()
            {
                if runtime_options.progress_enabled {
                    emit_progress(
                        "sketch_save",
                        &format!(
                            "event=pack_contigs\tcontigs_done={contigs_done}\tcontig_count={}\treference_minimizers_done={}",
                            contigs.len(),
                            reference_minimizer_count
                        ),
                        save_start,
                    );
                }
                check_memory_limit("sketch save while packing contigs", runtime_options)?;
            }
        }

        let contig_names: Vec<ReferenceContigName> =
            self.contig_names.clone().unwrap_or_else(|| {
                contigs
                    .iter()
                    .enumerate()
                    .map(|(contig_id, contig)| ReferenceContigName {
                        file_id: contig.file_id,
                        name: format!("contig_{contig_id}"),
                        segment_start: 0,
                        segment_end: 0,
                    })
                    .collect()
            });
        if contig_names.len() != contig_records.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "contig sidecar would have {} records, expected {}",
                    contig_names.len(),
                    contig_records.len()
                ),
            ));
        }
        write_contig_name_sidecar(
            path,
            &files,
            &contig_names,
            runtime_options.effective_worker_threads(),
        )?;

        let metadata: CachedReferenceMetadata = CachedReferenceMetadata {
            version: SKETCH_VERSION,
            k: kmer_size,
            w: window_size,
            key_mode: SKETCH_KEY_MODE.to_string(),
            fragment_length,
            min_fragment_length,
            split_n_run,
            dust_enabled: false,
            files,
            mphf,
            key_count: slot_keys.len(),
            hit_count: hit_payloads.len(),
            contig_count: contig_records.len(),
            reference_minimizer_count,
        };
        let metadata_bytes: Vec<u8> = serde_json::to_vec(&metadata).map_err(|err| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("failed to encode reference sketch metadata: {err}"),
            )
        })?;
        if runtime_options.progress_enabled {
            emit_progress(
                "sketch_save",
                &format!(
                    "event=metadata_encoded\tmetadata_bytes={}",
                    metadata_bytes.len()
                ),
                save_start,
            );
        }
        check_memory_limit("sketch save after metadata encode", runtime_options)?;
        let metadata_end: usize = SKETCH_MAGIC
            .len()
            .checked_add(size_of::<u64>())
            .and_then(|offset| offset.checked_add(metadata_bytes.len()))
            .ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "metadata length overflow")
            })?;
        let slot_keys_offset: usize = align_up(metadata_end, 8);
        let hit_offsets_offset: usize = align_up(
            checked_section_end(slot_keys_offset, slot_keys.len(), size_of::<MinimizerKey>())?,
            align_of::<u64>(),
        );
        let hit_counts_offset: usize =
            checked_section_end(hit_offsets_offset, hit_offsets.len(), size_of::<u64>())?;
        let hit_payloads_offset: usize = align_up(
            checked_section_end(hit_counts_offset, hit_counts.len(), size_of::<u32>())?,
            align_of::<SeedHit>(),
        );
        let contig_records_offset: usize = align_up(
            checked_section_end(
                hit_payloads_offset,
                hit_payloads.len(),
                size_of::<SeedHit>(),
            )?,
            align_of::<ContigRecord>(),
        );
        let reference_minimizers_offset: usize = align_up(
            checked_section_end(
                contig_records_offset,
                contig_records.len(),
                size_of::<ContigRecord>(),
            )?,
            align_of::<ReferenceMinimizer>(),
        );

        let mut sketch_output: SketchOutput = SketchOutput::create(
            path,
            tmp_dir,
            false,
            runtime_options.effective_worker_threads(),
        )?;
        let mut writer: &mut BufWriter<fs::File> = sketch_output.writer_mut()?;
        writer.write_all(SKETCH_MAGIC)?;
        writer.write_all(&(metadata_bytes.len() as u64).to_le_bytes())?;
        writer.write_all(&metadata_bytes)?;
        write_padding(&mut writer, slot_keys_offset - metadata_end)?;
        writer.write_all(slice_as_bytes(&slot_keys))?;
        write_padding(
            &mut writer,
            hit_offsets_offset
                - checked_section_end(
                    slot_keys_offset,
                    slot_keys.len(),
                    size_of::<MinimizerKey>(),
                )?,
        )?;
        writer.write_all(slice_as_bytes(&hit_offsets))?;
        writer.write_all(slice_as_bytes(&hit_counts))?;
        write_padding(
            &mut writer,
            hit_payloads_offset
                - checked_section_end(hit_counts_offset, hit_counts.len(), size_of::<u32>())?,
        )?;
        writer.write_all(slice_as_bytes(&hit_payloads))?;
        write_padding(
            &mut writer,
            contig_records_offset
                - checked_section_end(
                    hit_payloads_offset,
                    hit_payloads.len(),
                    size_of::<SeedHit>(),
                )?,
        )?;
        writer.write_all(slice_as_bytes(&contig_records))?;
        write_padding(
            &mut writer,
            reference_minimizers_offset
                - checked_section_end(
                    contig_records_offset,
                    contig_records.len(),
                    size_of::<ContigRecord>(),
                )?,
        )?;
        reference_minimizer_writer.flush()?;
        drop(reference_minimizer_writer);
        let expected_reference_minimizer_bytes: u64 =
            (reference_minimizer_count * size_of::<ReferenceMinimizer>()) as u64;
        let actual_reference_minimizer_bytes: u64 =
            fs::metadata(&reference_minimizer_scratch.path)?.len();
        if actual_reference_minimizer_bytes != expected_reference_minimizer_bytes {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "reference minimizer scratch file has {actual_reference_minimizer_bytes} bytes, expected {expected_reference_minimizer_bytes}"
                ),
            ));
        }
        let mut reference_minimizer_reader: fs::File =
            fs::File::open(&reference_minimizer_scratch.path)?;
        io::copy(&mut reference_minimizer_reader, &mut writer)?;
        let output_file_bytes: u64 = sketch_output.finish()?;
        if runtime_options.progress_enabled {
            emit_progress(
                "sketch_save",
                &format!(
                    "event=complete\tpath={}\tfile_bytes={}",
                    path.display(),
                    output_file_bytes
                ),
                save_start,
            );
        }
        check_memory_limit("sketch save complete", runtime_options)
    }

    /// Load a previously saved sketch cache and mmap its hit arrays.
    fn load(
        path: &Path,
        kmer_size: usize,
        window_size: usize,
        fragment_length: u32,
        min_fragment_length: u32,
        split_n_run: usize,
        load_contig_names: bool,
        tmp_dir: Option<&Path>,
        runtime_options: RuntimeOptions,
    ) -> io::Result<Self> {
        let load_start: Instant = Instant::now();
        if runtime_options.progress_enabled {
            emit_progress(
                "sketch_load",
                &format!("event=start\tpath={}", path.display()),
                load_start,
            );
        }
        check_memory_limit("sketch load start", runtime_options)?;

        let decompressed_sketch: Option<ScratchFile> = if is_gzip_path(path) {
            if runtime_options.progress_enabled {
                emit_progress(
                    "sketch_load",
                    &format!("event=decompress_start\tpath={}", path.display()),
                    load_start,
                );
            }
            Some(decompress_to_scratch(path, tmp_dir, "decompressed-sketch")?)
        } else {
            None
        };
        let mmap_path: &Path = decompressed_sketch
            .as_ref()
            .map(|scratch| scratch.path.as_path())
            .unwrap_or(path);

        let mut file: fs::File = fs::File::open(mmap_path)?;
        let mut magic: [u8; 8] = [0u8; 8];
        file.read_exact(&mut magic)?;
        if &magic != SKETCH_MAGIC {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "reference sketch cache {} has an invalid magic header",
                    path.display()
                ),
            ));
        }

        let mut metadata_len: [u8; 8] = [0u8; 8];
        file.read_exact(&mut metadata_len)?;
        let metadata_len: usize = u64::from_le_bytes(metadata_len) as usize;
        let metadata_start: usize = SKETCH_MAGIC.len() + size_of::<u64>();
        let metadata_end: usize = metadata_start.checked_add(metadata_len).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "metadata length overflow")
        })?;

        let mmap: Arc<MmapFile> = Arc::new(MmapFile::open(mmap_path)?);
        let bytes: &[u8] = mmap.as_slice();
        if metadata_end > bytes.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "reference sketch cache metadata extends past end of file",
            ));
        }

        let cached: CachedReferenceMetadata =
            serde_json::from_slice(&bytes[metadata_start..metadata_end]).map_err(|err| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("failed to read reference sketch cache: {err}"),
                )
            })?;
        if runtime_options.progress_enabled {
            emit_progress(
                "sketch_load",
                &format!(
                    "event=metadata_loaded\tfiles={}\tcontigs={}\tkey_count={}\thit_count={}\treference_minimizers={}",
                    cached.files.len(),
                    cached.contig_count,
                    cached.key_count,
                    cached.hit_count,
                    cached.reference_minimizer_count
                ),
                load_start,
            );
        }
        check_memory_limit("sketch load after metadata", runtime_options)?;

        if cached.dust_enabled {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "reference sketch cache is incompatible: it was built with the removed --dust filter",
            ));
        }

        if cached.version != SKETCH_VERSION
            || cached.k != kmer_size
            || cached.w != window_size
            || cached.key_mode != SKETCH_KEY_MODE
            || cached.fragment_length != fragment_length
            || cached.min_fragment_length != min_fragment_length
            || cached.split_n_run != split_n_run
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "reference sketch cache is incompatible: version={} k={} w={} key_mode={} fragment_length={} min_fragment_length={} split_n_run={}",
                    cached.version,
                    cached.k,
                    cached.w,
                    cached.key_mode,
                    cached.fragment_length,
                    cached.min_fragment_length,
                    cached.split_n_run
                ),
            ));
        }

        let slot_keys_offset: usize = align_up(metadata_end, 8);
        let hit_offsets_offset: usize = align_up(
            checked_section_end(
                slot_keys_offset,
                cached.key_count,
                size_of::<MinimizerKey>(),
            )?,
            align_of::<u64>(),
        );
        let hit_counts_offset: usize =
            checked_section_end(hit_offsets_offset, cached.key_count, size_of::<u64>())?;
        let hit_payloads_offset: usize = align_up(
            checked_section_end(hit_counts_offset, cached.key_count, size_of::<u32>())?,
            align_of::<SeedHit>(),
        );
        let contig_records_offset: usize = align_up(
            checked_section_end(hit_payloads_offset, cached.hit_count, size_of::<SeedHit>())?,
            align_of::<ContigRecord>(),
        );
        let reference_minimizers_offset: usize = align_up(
            checked_section_end(
                contig_records_offset,
                cached.contig_count,
                size_of::<ContigRecord>(),
            )?,
            align_of::<ReferenceMinimizer>(),
        );
        let file_end: usize = checked_section_end(
            reference_minimizers_offset,
            cached.reference_minimizer_count,
            size_of::<ReferenceMinimizer>(),
        )?;
        if file_end > bytes.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "reference sketch cache array sections extend past end of file",
            ));
        }

        if runtime_options.progress_enabled {
            emit_progress(
                "sketch_load",
                &format!(
                    "event=complete\tmmap_bytes={}\tfiles={}\tcontigs={}\tunique_minimizers={}",
                    bytes.len(),
                    cached.files.len(),
                    cached.contig_count,
                    cached.key_count
                ),
                load_start,
            );
        }
        check_memory_limit("sketch load complete", runtime_options)?;

        let contig_names: Option<Vec<ReferenceContigName>> = if load_contig_names {
            Some(load_contig_name_sidecar(path, cached.contig_count)?)
        } else {
            None
        };

        Ok(Self {
            files: cached.files,
            contigs: ReferenceContigs::Mmap(MmapReferenceContigs {
                mmap: Arc::clone(&mmap),
                contig_count: cached.contig_count,
                reference_minimizer_count: cached.reference_minimizer_count,
                contig_records_offset,
                reference_minimizers_offset,
            }),
            contig_names,
            index: ReferenceIndex::Mphf(MmapReferenceIndex {
                mphf: cached.mphf,
                mmap,
                key_count: cached.key_count,
                hit_count: cached.hit_count,
                slot_keys_offset,
                hit_offsets_offset,
                hit_counts_offset,
                hit_payloads_offset,
            }),
        })
    }

    /// Collect and merge seed-hit candidate intervals for one query fragment.
    fn find_candidate_regions(
        &self,
        query_minimizers: &[MinimizerKey],
        fragment_length: u32,
        minimum_shared_minimizers: usize,
        frequency_threshold: usize,
        seed_hits: &mut Vec<SeedHit>,
        candidate_regions: &mut Vec<ReferenceCandidateRegion>,
        #[cfg(debug_assertions)] mut mapping_metrics: Option<&mut MappingMetrics>,
    ) {
        seed_hits.clear();
        candidate_regions.clear();

        for minimizer in query_minimizers {
            let hits: Option<&[SeedHit]> = self.index.get(minimizer);
            #[cfg(debug_assertions)]
            if let Some(metrics) = mapping_metrics.as_deref_mut() {
                metrics.record_seed_lookup(hits.map(<[SeedHit]>::len), frequency_threshold);
            }

            if let Some(hits) = hits {
                if hits.len() < frequency_threshold {
                    seed_hits.extend_from_slice(hits);
                }
            }
        }

        seed_hits.sort_unstable_by_key(|hit| (hit.reference_contig_id, hit.position));

        let minimum_shared_minimizers: usize = minimum_shared_minimizers.max(1);

        for i in 0..seed_hits.len() {
            let Some(j) = i.checked_add(minimum_shared_minimizers - 1) else {
                break;
            };

            if j >= seed_hits.len() {
                break;
            }

            let first: SeedHit = seed_hits[i];
            let last: SeedHit = seed_hits[j];

            if first.reference_contig_id != last.reference_contig_id {
                continue;
            }

            if last.position.saturating_sub(first.position) >= fragment_length {
                continue;
            }

            let candidate_region: ReferenceCandidateRegion = ReferenceCandidateRegion {
                reference_contig_id: first.reference_contig_id as usize,
                start_position: last.position.saturating_sub(fragment_length - 1),
                end_position: first.position,
            };

            if let Some(previous) = candidate_regions.last_mut() {
                if previous.reference_contig_id == candidate_region.reference_contig_id
                    && previous.end_position >= candidate_region.start_position
                {
                    previous.end_position =
                        previous.end_position.max(candidate_region.end_position);
                    continue;
                }
            }

            candidate_regions.push(candidate_region);
        }
    }

    /// Score one candidate region with FastANI-style sliding-window minimizer overlap.
    #[cfg_attr(not(debug_assertions), allow(unused_mut, unused_variables))]
    fn score_candidate_region(
        &self,
        query_fragment: &QueryFragment,
        candidate_region: ReferenceCandidateRegion,
        kmer_size: usize,
        window_size: usize,
        min_identity: f64,
        mash_confidence: f64,
        counter: &mut SlidingSketchCounter,
        mut mapping_metrics: Option<&mut MappingMetrics>,
    ) -> Option<MappingResult> {
        #[cfg(not(debug_assertions))]
        let _ = mapping_metrics;

        let reference_file_id: usize =
            self.contigs.file_id(candidate_region.reference_contig_id)?;
        let minimizers: &[ReferenceMinimizer] = self
            .contigs
            .minimizers(candidate_region.reference_contig_id)?;

        if minimizers.is_empty() || query_fragment.minimizers.is_empty() {
            return None;
        }

        let count_minimizer_windows: u32 = query_fragment.length.saturating_sub(
            (window_size as u32).saturating_sub(1) + (kmer_size as u32).saturating_sub(1),
        );
        let first_start: usize =
            lower_bound_minimizer_position(minimizers, candidate_region.start_position);

        if first_start >= minimizers.len() {
            return None;
        }

        let last_end: usize = lower_bound_minimizer_position(
            minimizers,
            candidate_region
                .end_position
                .saturating_add(query_fragment.length),
        );
        #[cfg(debug_assertions)]
        {
            if let Some(metrics) = mapping_metrics.as_deref_mut() {
                metrics.reference_minimizers_scanned += last_end.saturating_sub(first_start);
            }
        }
        counter.prepare(
            &query_fragment.minimizers,
            &minimizers[first_start..last_end],
        );

        let mut best_shared: usize = 0usize;
        let mut best_reference_minimizer_count: usize = 0usize;
        let mut best_start: u32 = minimizers[first_start].position;
        let mut last_best_start: u32 = best_start;
        let mut window_end: usize = first_start;
        #[cfg(debug_assertions)]
        let mut scoring_window_steps: usize = 0usize;

        for start_idx in first_start..last_end {
            let start_position: u32 = minimizers[start_idx].position;

            if start_position > candidate_region.end_position {
                break;
            }
            #[cfg(debug_assertions)]
            {
                scoring_window_steps += 1;
            }

            if start_idx != first_start {
                counter.remove(minimizers[start_idx - 1].hash);
            }

            let end_position: u32 = start_position.saturating_add(count_minimizer_windows);

            while window_end < last_end && minimizers[window_end].position < end_position {
                counter.insert(minimizers[window_end].hash);
                window_end += 1;
            }

            let shared: usize = counter.shared_count();

            if shared > best_shared {
                best_shared = shared;
                best_reference_minimizer_count = counter.reference_minimizer_count();
                best_start = start_position;
                last_best_start = start_position;
            } else if shared == best_shared {
                best_reference_minimizer_count = counter.reference_minimizer_count();
                last_best_start = start_position;
            }
        }
        #[cfg(debug_assertions)]
        {
            if let Some(metrics) = mapping_metrics.as_deref_mut() {
                metrics.scoring_window_steps += scoring_window_steps;
            }
        }

        if best_shared == 0 {
            return None;
        }

        let sketch_size: usize = query_fragment.minimizers.len();
        let jaccard: f64 = best_shared as f64 / sketch_size as f64;
        let distance: f64 = fastani_mash_distance(jaccard, kmer_size);
        let lower_distance: f64 =
            mash_distance_lower_bound(distance, sketch_size, kmer_size, mash_confidence);
        let identity: f64 = 100.0 * (1.0 - distance);
        let upper_identity: f64 = 100.0 * (1.0 - lower_distance);

        if upper_identity < min_identity {
            return None;
        }

        Some(MappingResult {
            reference_file_id,
            reference_contig_id: candidate_region.reference_contig_id,
            query_fragment_id: query_fragment.id,
            query_fragment_length: query_fragment.length,
            reference_start: (best_start + last_best_start) / 2,
            identity,
            query_minimizer_count: sketch_size,
            reference_minimizer_count: best_reference_minimizer_count,
            shared_minimizers: best_shared,
            union_minimizers: sketch_size,
            jaccard,
        })
    }
}

impl SketchDatabase {
    fn collect_or_load(
        reference_paths: &[String],
        kmer_size: usize,
        window_size: usize,
        fragment_length: u32,
        min_fragment_length: u32,
        split_n_run: usize,
        sketch_prefix: Option<&Path>,
        tmp_dir: Option<&Path>,
        bgzip: bool,
        shard_size: usize,
        shard_minimizers: usize,
        index_build_mode: IndexBuildMode,
        threads: usize,
        load_contig_names: bool,
        runtime_options: RuntimeOptions,
    ) -> io::Result<Self> {
        let Some(prefix) = sketch_prefix else {
            return Ok(Self::Single(ReferenceSketch::collect(
                reference_paths,
                kmer_size,
                window_size,
                fragment_length,
                min_fragment_length,
                split_n_run,
                runtime_options,
            )?));
        };

        validate_shard_size(shard_size)?;
        validate_shard_minimizers(shard_minimizers)?;

        let manifest_path: PathBuf = manifest_path(prefix);
        if manifest_path.exists() {
            let manifest: ShardManifest = Self::load_manifest(
                prefix,
                kmer_size,
                window_size,
                fragment_length,
                min_fragment_length,
                split_n_run,
            )?;
            return Ok(Self::Sharded {
                prefix: prefix.to_path_buf(),
                manifest,
            });
        }

        if let Some(legacy_path) = legacy_sketch_path(prefix) {
            return Ok(Self::Single(ReferenceSketch::load(
                &legacy_path,
                kmer_size,
                window_size,
                fragment_length,
                min_fragment_length,
                split_n_run,
                load_contig_names,
                tmp_dir,
                runtime_options,
            )?));
        }

        let manifest: ShardManifest = Self::build_sharded(
            reference_paths,
            kmer_size,
            window_size,
            fragment_length,
            min_fragment_length,
            split_n_run,
            prefix,
            tmp_dir,
            bgzip,
            shard_size,
            shard_minimizers,
            index_build_mode,
            threads,
            runtime_options,
        )?;

        Ok(Self::Sharded {
            prefix: prefix.to_path_buf(),
            manifest,
        })
    }

    fn build_sharded(
        reference_paths: &[String],
        kmer_size: usize,
        window_size: usize,
        fragment_length: u32,
        min_fragment_length: u32,
        split_n_run: usize,
        prefix: &Path,
        tmp_dir: Option<&Path>,
        bgzip: bool,
        shard_size: usize,
        shard_minimizers: usize,
        index_build_mode: IndexBuildMode,
        threads: usize,
        runtime_options: RuntimeOptions,
    ) -> io::Result<ShardManifest> {
        validate_shard_size(shard_size)?;
        validate_shard_minimizers(shard_minimizers)?;
        if let Some(parent) = prefix
            .parent()
            .filter(|parent| !parent.as_os_str().is_empty())
        {
            fs::create_dir_all(parent)?;
        }

        let build_start: Instant = Instant::now();
        if runtime_options.progress_enabled {
            emit_progress(
                "database_build",
                &format!(
                    "event=start\tmode=sharded\tprefix={}\treferences={}\tshard_size={shard_size}\tshard_minimizers={shard_minimizers}\tthreads={}",
                    prefix.display(),
                    reference_paths.len(),
                    threads
                ),
                build_start,
            );
        }

        let shard_plans: Vec<ShardPlan> = plan_shards_by_minimizers(
            reference_paths,
            kmer_size,
            window_size,
            split_n_run,
            shard_size,
            shard_minimizers,
            threads,
            runtime_options,
        )?;
        let build_parallelism: usize =
            database_build_parallelism(threads, &shard_plans, runtime_options.max_memory_bytes);

        if runtime_options.progress_enabled {
            emit_progress(
                "database_build",
                &format!(
                    "event=shards_planned\tshards={}\tbuild_parallelism={build_parallelism}\tthreads={threads}\tmax_memory_gb={}\tshard_minimizers={shard_minimizers}",
                    shard_plans.len(),
                    runtime_options
                        .max_memory_bytes
                        .map(|bytes| format!("{:.3}", bytes as f64 / (1024.0 * 1024.0 * 1024.0)))
                        .unwrap_or_else(|| "unset".to_string())
                ),
                build_start,
            );
        }

        let completed_shards: Arc<AtomicUsize> = Arc::new(AtomicUsize::new(0));
        let shard_worker_threads: usize = threads.max(1).div_ceil(build_parallelism.max(1));
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(build_parallelism)
            .build()
            .map_err(|err| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("failed to initialize database build thread pool: {err}"),
                )
            })?;

        let mut shard_results: Vec<ShardBuildResult> = pool.install(|| {
            shard_plans
                .par_iter()
                .enumerate()
                .map(|(shard_offset, shard_plan)| {
                    let shard_index: usize = shard_offset + 1;
                    let first_reference: usize = shard_plan.first_reference;
                    let shard_end: usize = first_reference
                        .checked_add(shard_plan.reference_count)
                        .ok_or_else(|| {
                            io::Error::new(
                                io::ErrorKind::InvalidData,
                                "shard reference range overflow",
                            )
                        })?;
                    let reference_chunk: &[String] = &reference_paths[first_reference..shard_end];
                    let shard_path: PathBuf = shard_path(prefix, shard_index, bgzip);

                    if runtime_options.progress_enabled {
                        let effective_index_build_mode: IndexBuildMode =
                            effective_index_build_mode(
                                index_build_mode,
                                shard_plan.estimated_minimizers,
                            );
                        emit_progress(
                            "database_build",
                            &format!(
                                "event=shard_start\tshard={shard_index}\tfirst_reference={first_reference}\treference_count={}\testimated_minimizers={}\testimated_memory_mib={:.3}\tindex_build_mode={}\trequested_index_build_mode={}\tpath={}",
                                reference_chunk.len(),
                                shard_plan.estimated_minimizers,
                                memory_mib(estimate_partitioned_shard_memory_bytes(
                                    shard_plan.estimated_minimizers
                                )),
                                effective_index_build_mode.name(),
                                index_build_mode.name(),
                                shard_path.display()
                            ),
                            build_start,
                        );
                    }

                    let stats: SketchBuildStats = ReferenceSketch::collect_and_save_streaming(
                        reference_chunk,
                        kmer_size,
                        window_size,
                        fragment_length,
                        min_fragment_length,
                        split_n_run,
                        &shard_path,
                        tmp_dir,
                        bgzip,
                        shard_plan.estimated_minimizers,
                        index_build_mode,
                        runtime_options.with_worker_threads(shard_worker_threads),
                    )?;

                    let shards_done: usize =
                        completed_shards.fetch_add(1, AtomicOrdering::Relaxed) + 1;
                    if runtime_options.progress_enabled {
                        emit_progress(
                            "database_build",
                            &format!(
                                "event=shard_complete\tshard={shard_index}\tshards_done={shards_done}\tshards_total={}\treferences_done={}\treferences_total={}\tcontigs={}\treference_minimizers={}",
                                shard_plans.len(),
                                first_reference + stats.reference_count,
                                reference_paths.len(),
                                stats.reference_contig_count,
                                stats.reference_minimizer_count
                            ),
                            build_start,
                        );
                    }
                    check_memory_limit("after sharded sketch build shard", runtime_options)?;

                    Ok(ShardBuildResult {
                        entry: ShardManifestEntry {
                            shard_index,
                            filename: shard_filename(prefix, shard_index, bgzip),
                            first_reference,
                            reference_count: stats.reference_count,
                            reference_contigs: stats.reference_contig_count,
                            mapped_reference_length: stats.mapped_reference_length,
                            reference_minimizers: stats.reference_minimizer_count,
                            unique_minimizers: stats.unique_minimizer_count,
                        },
                    })
                })
                .collect::<io::Result<Vec<_>>>()
        })?;

        shard_results.sort_by_key(|result| result.entry.shard_index);
        let shards: Vec<ShardManifestEntry> = shard_results
            .into_iter()
            .map(|result| result.entry)
            .collect();
        let total_reference_contigs: usize =
            shards.iter().map(|shard| shard.reference_contigs).sum();
        let total_mapped_reference_length: u64 = shards
            .iter()
            .map(|shard| shard.mapped_reference_length)
            .sum();
        let total_reference_minimizers: usize =
            shards.iter().map(|shard| shard.reference_minimizers).sum();
        let total_shard_unique_minimizers: usize =
            shards.iter().map(|shard| shard.unique_minimizers).sum();

        let manifest: ShardManifest = ShardManifest {
            sketch_format_version: SKETCH_VERSION,
            database_schema_version: SKETCH_DATABASE_SCHEMA_VERSION,
            k: kmer_size,
            w: window_size,
            key_mode: SKETCH_KEY_MODE.to_string(),
            fragment_length,
            min_fragment_length,
            split_n_run,
            dust_enabled: false,
            shard_size,
            shard_minimizers,
            total_references: reference_paths.len(),
            total_reference_contigs,
            total_mapped_reference_length,
            total_reference_minimizers,
            total_shard_unique_minimizers,
            build_unix_seconds: unix_timestamp_seconds()?,
            build_args: env::args().collect(),
            reference_list_checksum: reference_list_checksum(reference_paths),
            shards,
        };

        Self::write_manifest(prefix, &manifest)?;
        if runtime_options.progress_enabled {
            emit_progress(
                "database_build",
                &format!(
                    "event=complete\tmanifest={}\tshards={}\treferences={}\tcontigs={}\treference_minimizers={}",
                    manifest_path(prefix).display(),
                    manifest.shards.len(),
                    manifest.total_references,
                    manifest.total_reference_contigs,
                    manifest.total_reference_minimizers
                ),
                build_start,
            );
        }

        Ok(manifest)
    }

    fn write_manifest(prefix: &Path, manifest: &ShardManifest) -> io::Result<()> {
        let path: PathBuf = manifest_path(prefix);
        if let Some(parent) = path
            .parent()
            .filter(|parent| !parent.as_os_str().is_empty())
        {
            fs::create_dir_all(parent)?;
        }

        let manifest_bytes: Vec<u8> = serde_json::to_vec_pretty(manifest).map_err(|err| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("failed to encode sharded sketch manifest: {err}"),
            )
        })?;
        fs::write(path, manifest_bytes)
    }

    fn load_manifest(
        prefix: &Path,
        kmer_size: usize,
        window_size: usize,
        fragment_length: u32,
        min_fragment_length: u32,
        split_n_run: usize,
    ) -> io::Result<ShardManifest> {
        let path: PathBuf = manifest_path(prefix);
        let manifest_bytes: Vec<u8> = fs::read(&path)?;
        let manifest: ShardManifest = serde_json::from_slice(&manifest_bytes).map_err(|err| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "failed to read sharded sketch manifest {}: {err}",
                    path.display()
                ),
            )
        })?;

        if let Some(error) = shard_manifest_compatibility_error(
            &manifest,
            kmer_size,
            window_size,
            fragment_length,
            min_fragment_length,
            split_n_run,
        ) {
            return Err(io::Error::new(io::ErrorKind::InvalidData, error));
        }

        let reference_count_sum: usize = manifest
            .shards
            .iter()
            .map(|shard| shard.reference_count)
            .sum();
        if reference_count_sum != manifest.total_references {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "sharded sketch manifest reference count mismatch: shards={reference_count_sum} total={}",
                    manifest.total_references
                ),
            ));
        }

        for shard in &manifest.shards {
            let path: PathBuf = shard_entry_path(prefix, shard);
            if !path.exists() {
                return Err(io::Error::new(
                    io::ErrorKind::NotFound,
                    format!(
                        "sharded sketch manifest references missing shard: {}",
                        path.display()
                    ),
                ));
            }
        }

        Ok(manifest)
    }

    fn reference_count(&self) -> usize {
        match self {
            Self::Single(sketch) => sketch.files.len(),
            Self::Sharded { manifest, .. } => manifest.total_references,
        }
    }

    fn contig_count(&self) -> usize {
        match self {
            Self::Single(sketch) => sketch.contigs.len(),
            Self::Sharded { manifest, .. } => manifest.total_reference_contigs,
        }
    }

    fn unique_minimizer_count(&self) -> usize {
        match self {
            Self::Single(sketch) => sketch.index.len(),
            Self::Sharded { manifest, .. } => manifest.total_shard_unique_minimizers,
        }
    }

    fn mode_name(
        &self,
        sketch_was_requested: bool,
        existing_database_loaded: bool,
    ) -> &'static str {
        match self {
            Self::Single(_) if !sketch_was_requested => "build-fasta",
            Self::Single(_) if existing_database_loaded => "load-sketch",
            Self::Single(_) => "build-fasta",
            Self::Sharded { .. } if existing_database_loaded => "load-sharded-sketch",
            Self::Sharded { .. } => "build-sharded-sketch",
        }
    }
}

impl QueryFile {
    /// Read a query FASTA file and split each contig into query fragments.
    fn collect(
        reader: &mut fasta::io::Reader<impl io::BufRead>,
        kmer_size: usize,
        window_size: usize,
        minmer_count: Option<usize>,
        fragment_length: u32,
        fragment_stride: u32,
        min_fragment_length: u32,
        split_n_run: usize,
    ) -> io::Result<Self> {
        let mut fragments: Vec<QueryFragment> = Vec::new();
        let mut contig_names: Vec<String> = Vec::new();
        let mut mapped_length: u64 = 0u64;
        let fragment_length: usize = fragment_length as usize;
        let fragment_stride: usize = fragment_stride as usize;
        let min_fragment_length: usize = min_fragment_length as usize;

        for result in reader.records() {
            let record = result.map_err(|err| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("failed to read FASTA record from query: {err}"),
                )
            })?;
            let contig_id: usize = contig_names.len();
            contig_names.push(String::from_utf8_lossy(record.name()).into_owned());
            let sequence: &fasta::record::Sequence = record.sequence();
            let sequence_bytes: &[u8] = sequence.as_ref();

            for segment_range in split_sequence_ranges(sequence_bytes, split_n_run) {
                let segment_start: usize = segment_range.start;
                let segment_sequence: &[u8] = &sequence_bytes[segment_range];
                let fragment_ranges: Vec<std::ops::Range<usize>> = query_fragment_ranges(
                    segment_sequence.len(),
                    fragment_length,
                    fragment_stride,
                    min_fragment_length,
                );
                fragments.reserve(fragment_ranges.len());

                for fragment_range in fragment_ranges {
                    let fragment_sketch: QueryFragmentSketch = query_fragment_sketch(
                        &segment_sequence[fragment_range.clone()],
                        kmer_size,
                        window_size,
                        minmer_count,
                    );
                    if fragment_sketch.minimizers.is_empty() {
                        continue;
                    }

                    let fragment_length: u32 = fragment_range.len() as u32;
                    let query_start: u32 = (segment_start + fragment_range.start) as u32;
                    let query_end: u32 = (segment_start + fragment_range.end) as u32;
                    mapped_length += u64::from(fragment_length);
                    fragments.push(QueryFragment {
                        id: fragments.len(),
                        contig_id,
                        start: query_start,
                        end: query_end,
                        length: fragment_length,
                        minimizers: fragment_sketch.minimizers,
                        seed_minimizers: fragment_sketch.seed_minimizers,
                    });
                }
            }
        }

        if fragments.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "ERROR: Input has no usable fragments",
            ));
        }

        Ok(Self {
            fragments,
            contig_names,
            mapped_length,
        })
    }

    fn mapped_length(&self) -> u64 {
        self.mapped_length
    }

    fn total_minimizers(&self) -> usize {
        self.fragments
            .iter()
            .map(|fragment| fragment.minimizers.len())
            .sum()
    }

    fn total_seed_minimizers(&self) -> usize {
        self.fragments
            .iter()
            .map(|fragment| fragment.seed_minimizers.len())
            .sum()
    }

    #[cfg(debug_assertions)]
    fn memory_estimate(&self) -> QueryMemoryEstimate {
        QueryMemoryEstimate {
            fragment_struct_bytes: self.fragments.capacity() * size_of::<QueryFragment>(),
            query_minimizer_vec_bytes: self
                .fragments
                .iter()
                .map(|fragment| fragment.minimizers.capacity() * size_of::<MinimizerKey>())
                .sum(),
            seed_minimizer_vec_bytes: self
                .fragments
                .iter()
                .map(|fragment| fragment.seed_minimizers.capacity() * size_of::<MinimizerKey>())
                .sum(),
        }
    }
}

/// Map all query fragments to the reference sketch on one thread.
fn map_query_to_reference(
    reference_sketch: &ReferenceSketch,
    query_file: &QueryFile,
    kmer_size: usize,
    window_size: usize,
    min_identity: f64,
    mash_confidence: f64,
    frequency_threshold: usize,
    collect_metrics: bool,
) -> MappingOutput {
    let mut mapping_results: Vec<MappingResult> = Vec::new();
    let mut scratch: MappingScratch = MappingScratch::default();
    let mut mapping_metrics: MappingMetrics = MappingMetrics::default();

    for query_fragment in &query_file.fragments {
        map_query_fragment_into(
            reference_sketch,
            query_fragment,
            kmer_size,
            window_size,
            min_identity,
            mash_confidence,
            frequency_threshold,
            &mut scratch,
            &mut mapping_results,
            &mut mapping_metrics,
            collect_metrics,
        );
    }

    MappingOutput {
        results: mapping_results,
        #[cfg(debug_assertions)]
        metrics: mapping_metrics,
    }
}

/// Map all query fragments to the reference sketch using a local Rayon thread pool.
fn map_query_to_reference_parallel(
    reference_sketch: &ReferenceSketch,
    query_file: &QueryFile,
    kmer_size: usize,
    window_size: usize,
    min_identity: f64,
    mash_confidence: f64,
    threads: usize,
    frequency_threshold: usize,
    collect_metrics: bool,
) -> io::Result<MappingOutput> {
    if threads <= 1 {
        return Ok(map_query_to_reference(
            reference_sketch,
            query_file,
            kmer_size,
            window_size,
            min_identity,
            mash_confidence,
            frequency_threshold,
            collect_metrics,
        ));
    }

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(|err| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("failed to initialize rayon thread pool: {err}"),
            )
        })?;

    Ok(pool.install(|| {
        query_file
            .fragments
            .par_iter()
            .fold(
                || {
                    (
                        MappingScratch::default(),
                        Vec::new(),
                        MappingMetrics::default(),
                    )
                },
                |(mut scratch, mut mapping_results, mut mapping_metrics), query_fragment| {
                    map_query_fragment_into(
                        reference_sketch,
                        query_fragment,
                        kmer_size,
                        window_size,
                        min_identity,
                        mash_confidence,
                        frequency_threshold,
                        &mut scratch,
                        &mut mapping_results,
                        &mut mapping_metrics,
                        collect_metrics,
                    );
                    (scratch, mapping_results, mapping_metrics)
                },
            )
            .map(|(_scratch, mapping_results, mapping_metrics)| {
                #[cfg(not(debug_assertions))]
                let _ = mapping_metrics;
                MappingOutput {
                    results: mapping_results,
                    #[cfg(debug_assertions)]
                    metrics: mapping_metrics,
                }
            })
            .reduce(
                || MappingOutput {
                    results: Vec::new(),
                    #[cfg(debug_assertions)]
                    metrics: MappingMetrics::default(),
                },
                |mut left, mut right| {
                    left.results.append(&mut right.results);
                    #[cfg(debug_assertions)]
                    left.metrics.merge(right.metrics);
                    left
                },
            )
    }))
}

#[cfg(debug_assertions)]
fn record_candidate_discovery_metrics(
    mapping_metrics: &mut MappingMetrics,
    seed_hit_count: usize,
    candidate_region_count: usize,
    elapsed: std::time::Duration,
) {
    mapping_metrics.candidate_discovery_calls += 1;
    mapping_metrics.seed_hits_collected += seed_hit_count;
    mapping_metrics.candidate_regions_found += candidate_region_count;
    mapping_metrics.candidate_discovery_elapsed += elapsed;
}

/// Map one query fragment, appending all reportAll-style surviving mappings.
fn map_query_fragment_into(
    reference_sketch: &ReferenceSketch,
    query_fragment: &QueryFragment,
    kmer_size: usize,
    window_size: usize,
    min_identity: f64,
    mash_confidence: f64,
    frequency_threshold: usize,
    scratch: &mut MappingScratch,
    mapping_results: &mut Vec<MappingResult>,
    mapping_metrics: &mut MappingMetrics,
    collect_metrics: bool,
) {
    #[cfg(not(debug_assertions))]
    {
        let _ = mapping_metrics;
        let _ = collect_metrics;
    }

    if query_fragment.seed_minimizers.is_empty() || query_fragment.minimizers.is_empty() {
        return;
    }

    let minimum_shared_seed_minimizers: usize = estimate_relaxed_minimum_shared_minimizers(
        query_fragment.seed_minimizers.len(),
        kmer_size,
        min_identity,
        mash_confidence,
    );
    #[cfg(debug_assertions)]
    let candidate_discovery_start: Option<Instant> = collect_metrics.then(Instant::now);
    reference_sketch.find_candidate_regions(
        &query_fragment.seed_minimizers,
        query_fragment.length,
        minimum_shared_seed_minimizers,
        frequency_threshold,
        &mut scratch.seed_hits,
        &mut scratch.candidate_regions,
        #[cfg(debug_assertions)]
        collect_metrics.then_some(&mut *mapping_metrics),
    );
    #[cfg(debug_assertions)]
    if let Some(start) = candidate_discovery_start {
        record_candidate_discovery_metrics(
            mapping_metrics,
            scratch.seed_hits.len(),
            scratch.candidate_regions.len(),
            start.elapsed(),
        );
    }

    scratch.fragment_mappings.clear();
    #[cfg(debug_assertions)]
    if collect_metrics {
        mapping_metrics.candidate_regions_scored += scratch.candidate_regions.len();
    }

    for &candidate_region in &scratch.candidate_regions {
        #[cfg(debug_assertions)]
        let scoring_start: Option<Instant> = collect_metrics.then(Instant::now);
        #[cfg(debug_assertions)]
        let mapping_metrics_ref: Option<&mut MappingMetrics> =
            collect_metrics.then_some(&mut *mapping_metrics);
        #[cfg(not(debug_assertions))]
        let mapping_metrics_ref: Option<&mut MappingMetrics> = None;
        let mapping: Option<MappingResult> = reference_sketch.score_candidate_region(
            query_fragment,
            candidate_region,
            kmer_size,
            window_size,
            min_identity,
            mash_confidence,
            &mut scratch.counter,
            mapping_metrics_ref,
        );
        #[cfg(debug_assertions)]
        if let Some(start) = scoring_start {
            mapping_metrics.scoring_elapsed += start.elapsed();
        }

        if let Some(mapping) = mapping {
            scratch.fragment_mappings.push(mapping);
        }
    }
    #[cfg(debug_assertions)]
    if collect_metrics {
        mapping_metrics.retained_mappings += scratch.fragment_mappings.len();
    }

    mapping_results.append(&mut scratch.fragment_mappings);
}

/// Collapse raw fragment mappings into final per-reference-file ANI summaries.
fn final_ani_computation(
    mut mapping_results: Vec<MappingResult>,
    reference_file_count: usize,
    fragment_length: u32,
    disable_reciprocal: bool,
) -> AniComputation {
    mapping_results.sort_by(compare_query_bucket);

    let mut query_best_mappings: Vec<MappingResult> = Vec::new();

    for mapping in mapping_results {
        if let Some(previous) = query_best_mappings.last_mut() {
            if previous.reference_file_id == mapping.reference_file_id
                && previous.query_fragment_id == mapping.query_fragment_id
            {
                *previous = mapping;
                continue;
            }
        }

        query_best_mappings.push(mapping);
    }

    let summary_mappings: Vec<MappingResult> = if disable_reciprocal {
        query_best_mappings
    } else {
        query_best_mappings
            .sort_by(|left, right| compare_refbin_bucket(left, right, fragment_length));

        let mut reciprocal_best_mappings: Vec<MappingResult> = Vec::new();

        for mapping in query_best_mappings {
            if let Some(previous) = reciprocal_best_mappings.last_mut() {
                if previous.reference_contig_id == mapping.reference_contig_id
                    && reference_position_bin(previous.reference_start, fragment_length)
                        == reference_position_bin(mapping.reference_start, fragment_length)
                {
                    *previous = mapping;
                    continue;
                }
            }

            reciprocal_best_mappings.push(mapping);
        }

        reciprocal_best_mappings
    };

    let mut summaries: Vec<AniSummary> = (0..reference_file_count)
        .map(|_| AniSummary::default())
        .collect::<Vec<_>>();
    let mut reciprocal_best_keys: HashSet<MappingResultKey> = HashSet::new();

    for mapping in summary_mappings {
        reciprocal_best_keys.insert(MappingResultKey::from_mapping(&mapping));
        let summary = &mut summaries[mapping.reference_file_id];
        summary.shared_fragments += 1;
        summary.shared_bases += u64::from(mapping.query_fragment_length);
        summary.weighted_identity_sum +=
            mapping.identity * f64::from(mapping.query_fragment_length);
    }

    AniComputation {
        summaries,
        reciprocal_best_keys,
    }
}

fn compare_query_bucket(left: &MappingResult, right: &MappingResult) -> Ordering {
    (
        left.reference_file_id,
        left.query_fragment_id,
        ordered_float(left.identity),
        left.reference_contig_id,
        left.reference_start,
    )
        .cmp(&(
            right.reference_file_id,
            right.query_fragment_id,
            ordered_float(right.identity),
            right.reference_contig_id,
            right.reference_start,
        ))
}

fn compare_refbin_bucket(
    left: &MappingResult,
    right: &MappingResult,
    fragment_length: u32,
) -> Ordering {
    (
        left.reference_contig_id,
        reference_position_bin(left.reference_start, fragment_length),
        ordered_float(left.identity),
    )
        .cmp(&(
            right.reference_contig_id,
            reference_position_bin(right.reference_start, fragment_length),
            ordered_float(right.identity),
        ))
}

fn ordered_float(value: f64) -> u64 {
    value.to_bits()
}

fn reference_position_bin(position: u32, fragment_length: u32) -> u32 {
    position / fragment_length.saturating_sub(20).max(1)
}

/// Parsed command-line arguments.
struct CliArgs {
    references: Vec<String>,
    queries: Vec<String>,
    sketch_path: Option<PathBuf>,
    tmp_dir: Option<PathBuf>,
    out_path: Option<PathBuf>,
    mapping_stats_path: Option<PathBuf>,
    bgzip: bool,
    verbose: bool,
    threads: usize,
    freq_threshold_percent: f64,
    minmer_count: Option<usize>,
    kmer_size: usize,
    window_size: usize,
    fragment_length: u32,
    fragment_stride: u32,
    min_fragment_length: u32,
    min_identity: f64,
    mash_confidence: f64,
    disable_reciprocal: bool,
    split_n_run: usize,
    max_memory_bytes: Option<u64>,
    shard_size: usize,
    shard_minimizers: usize,
    index_build_mode: IndexBuildMode,
}

fn usage() -> &'static str {
    "usage: fasterANI (--reference <reference.fa> | --reference-list <refs.txt>)... (--query <query.fa> | --query-list <queries.txt>)... [--sketch <prefix>] [--bgzip] [--kmer-size <n, default 16>] [--window-size <n, default 24>] [--fragment-length <bp, default 3000>] [--fragment-stride <bp, default fragment-length>] [--min-fraglen <bp, default fragment-length>] [--min-fragment-length <bp, default fragment-length>] [--min-identity <percent, default 80>] [--mash-confidence <0..1, default 0.9>] [--disable-reciprocal] [--shard-size <n, default 10000>] [--shard-minimizers <n, default memory-aware>] [--index-build-mode auto|hash|partitioned] [--tmp <dir>] [--out <output.tsv>] [--mapping-stats <output.tsv>] [--threads <n>] [--max-memory-gb <gb>] [--freq-threshold-percent <0..100>] [--minmer-count <n>] [--split-N <bp>] [--verbose]"
}

fn read_path_list(path: &str) -> io::Result<Vec<String>> {
    let contents: String = fs::read_to_string(path)?;
    Ok(contents
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(ToOwned::to_owned)
        .collect())
}

/// Parse command-line arguments.
fn parse_cli_args() -> io::Result<Option<CliArgs>> {
    let mut references: Vec<String> = Vec::new();
    let mut queries: Vec<String> = Vec::new();
    let mut sketch_path: Option<PathBuf> = None;
    let mut tmp_dir: Option<PathBuf> = None;
    let mut out_path: Option<PathBuf> = None;
    let mut mapping_stats_path: Option<PathBuf> = None;
    let mut bgzip: bool = false;
    let mut verbose: bool = false;
    let mut threads: usize = 1usize;
    let mut freq_threshold_percent: f64 = DEFAULT_FREQ_THRESHOLD_PERCENT;
    let mut minmer_count: Option<usize> = None;
    let mut kmer_size: usize = DEFAULT_KMER_SIZE;
    let mut window_size: usize = DEFAULT_WINDOW_SIZE;
    let mut fragment_length: u32 = DEFAULT_FRAGMENT_LENGTH;
    let mut fragment_stride: u32 = DEFAULT_FRAGMENT_STRIDE;
    let mut fragment_stride_was_set: bool = false;
    let mut min_fragment_length: u32 = DEFAULT_MIN_FRAGMENT_LENGTH;
    let mut min_fragment_length_was_set: bool = false;
    let mut min_identity: f64 = DEFAULT_MIN_PERCENT_IDENTITY;
    let mut mash_confidence: f64 = DEFAULT_MASH_CONFIDENCE;
    let mut disable_reciprocal: bool = false;
    let mut split_n_run: usize = DEFAULT_SPLIT_N_RUN;
    let mut max_memory_bytes: Option<u64> = None;
    let mut shard_size: usize = DEFAULT_SHARD_SIZE;
    let mut shard_minimizers: Option<usize> = None;
    let mut index_build_mode: IndexBuildMode = IndexBuildMode::Auto;
    let mut args: std::iter::Skip<std::env::Args> = env::args().skip(1);

    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--reference" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidInput, "--reference requires a path")
                })?;
                references.push(value);
            }
            "--reference-list" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--reference-list requires a path",
                    )
                })?;
                references.extend(read_path_list(&value)?);
            }
            "--query" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidInput, "--query requires a path")
                })?;
                queries.push(value);
            }
            "--query-list" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidInput, "--query-list requires a path")
                })?;
                queries.extend(read_path_list(&value)?);
            }
            "--sketch" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidInput, "--sketch requires a path")
                })?;
                sketch_path = Some(PathBuf::from(value));
            }
            "--bgzip" => {
                bgzip = true;
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
                validate_window_size(window_size)?;
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
                validate_fragment_length(fragment_length)?;
            }
            "--min-identity" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--min-identity requires a value",
                    )
                })?;
                min_identity = value.parse::<f64>().map_err(|err| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("invalid --min-identity value {value:?}: {err}"),
                    )
                })?;
                validate_min_identity(min_identity)?;
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
                validate_mash_confidence(mash_confidence)?;
            }
            "--disable-reciprocal" => {
                disable_reciprocal = true;
            }
            "--shard-size" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidInput, "--shard-size requires a value")
                })?;
                shard_size = value.parse::<usize>().map_err(|err| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("invalid --shard-size value {value:?}: {err}"),
                    )
                })?;
                validate_shard_size(shard_size)?;
            }
            "--shard-minimizers" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--shard-minimizers requires a value",
                    )
                })?;
                let parsed_shard_minimizers: usize = value.parse::<usize>().map_err(|err| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("invalid --shard-minimizers value {value:?}: {err}"),
                    )
                })?;
                validate_shard_minimizers(parsed_shard_minimizers)?;
                shard_minimizers = Some(parsed_shard_minimizers);
            }
            "--tmp" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidInput, "--tmp requires a directory")
                })?;
                tmp_dir = Some(PathBuf::from(value));
            }
            "--index-build-mode" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--index-build-mode requires a value",
                    )
                })?;
                index_build_mode = IndexBuildMode::parse(&value)?;
            }
            "--out" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidInput, "--out requires a path")
                })?;
                out_path = Some(PathBuf::from(value));
            }
            "--mapping-stats" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--mapping-stats requires a path",
                    )
                })?;
                mapping_stats_path = Some(PathBuf::from(value));
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
            }
            "--max-memory-gb" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--max-memory-gb requires a value",
                    )
                })?;
                max_memory_bytes = Some(parse_max_memory_gb(&value)?);
            }
            "--freq-threshold-percent" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--freq-threshold-percent requires a value",
                    )
                })?;
                freq_threshold_percent = value.parse::<f64>().map_err(|err| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("invalid --freq-threshold-percent value {value:?}: {err}"),
                    )
                })?;
                if !(0.0..=100.0).contains(&freq_threshold_percent) {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--freq-threshold-percent must be between 0 and 100",
                    ));
                }
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
                if min_fragment_length == 0 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("{arg} must be at least 1"),
                    ));
                }
            }
            "--split-N" => {
                let value = args.next().ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidInput, "--split-N requires a value")
                })?;
                split_n_run = value.parse::<usize>().map_err(|err| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("invalid --split-N value {value:?}: {err}"),
                    )
                })?;
            }
            "--verbose" => {
                verbose = true;
            }
            "--help" | "-h" => {
                eprintln!("{}", usage());
                return Ok(None);
            }
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("unknown argument {arg:?}\n{}", usage()),
                ));
            }
        }
    }

    if references.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("missing --reference\n{}", usage()),
        ));
    }

    if queries.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("missing --query\n{}", usage()),
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
    validate_min_identity(min_identity)?;
    validate_mash_confidence(mash_confidence)?;
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

    let shard_minimizers: usize = shard_minimizers
        .unwrap_or_else(|| default_shard_minimizers_for_runtime(threads, max_memory_bytes));
    validate_shard_minimizers(shard_minimizers)?;

    Ok(Some(CliArgs {
        references,
        queries,
        sketch_path,
        tmp_dir,
        out_path,
        mapping_stats_path,
        bgzip,
        verbose,
        threads,
        freq_threshold_percent,
        minmer_count,
        kmer_size,
        window_size,
        fragment_length,
        fragment_stride,
        min_fragment_length,
        min_identity,
        mash_confidence,
        disable_reciprocal,
        split_n_run,
        max_memory_bytes,
        shard_size,
        shard_minimizers,
        index_build_mode,
    }))
}

/// Return the current process peak RSS in kilobytes, or -1 when unavailable.
fn peak_rss_kb() -> i64 {
    let mut usage = std::mem::MaybeUninit::<libc::rusage>::uninit();
    let result = unsafe { libc::getrusage(libc::RUSAGE_SELF, usage.as_mut_ptr()) };
    if result == 0 {
        unsafe { usage.assume_init().ru_maxrss }
    } else {
        -1
    }
}

fn memory_mib(bytes: usize) -> f64 {
    bytes as f64 / 1024.0 / 1024.0
}

fn memory_gib(bytes: u64) -> f64 {
    bytes as f64 / 1024.0 / 1024.0 / 1024.0
}

fn current_rss_kb() -> i64 {
    let Ok(status) = fs::read_to_string("/proc/self/status") else {
        return peak_rss_kb();
    };

    for line in status.lines() {
        let Some(rest) = line.strip_prefix("VmRSS:") else {
            continue;
        };
        let Some(value) = rest.split_whitespace().next() else {
            continue;
        };
        if let Ok(kb) = value.parse::<i64>() {
            return kb;
        }
    }

    peak_rss_kb()
}

fn current_rss_bytes() -> Option<u64> {
    let rss_kb: i64 = current_rss_kb();
    if rss_kb <= 0 {
        None
    } else {
        Some(rss_kb as u64 * 1024)
    }
}

fn check_memory_limit(context: &str, options: RuntimeOptions) -> io::Result<()> {
    let Some(limit_bytes) = options.max_memory_bytes else {
        return Ok(());
    };
    let Some(rss_bytes) = current_rss_bytes() else {
        return Ok(());
    };

    if rss_bytes > limit_bytes {
        return Err(io::Error::new(
            io::ErrorKind::OutOfMemory,
            format!(
                "ERROR: {context} exceeded --max-memory-gb: rss_gb={:.3} limit_gb={:.3}",
                memory_gib(rss_bytes),
                memory_gib(limit_bytes)
            ),
        ));
    }

    Ok(())
}

fn emit_progress(stage: &str, message: &str, start: Instant) {
    let rss_kb: i64 = current_rss_kb();
    let rss_gib: f64 = if rss_kb > 0 {
        rss_kb as f64 / 1024.0 / 1024.0
    } else {
        f64::NAN
    };
    eprintln!(
        "PROGRESS\tstage={stage}\t{message}\trss_gib={rss_gib:.3}\telapsed_s={:.3}",
        start.elapsed().as_secs_f64()
    );
}

#[cfg(debug_assertions)]
fn reference_build_struct_bytes(
    reference_minimizers: usize,
    seed_hits: usize,
    unique_index_keys: usize,
) -> usize {
    reference_minimizers
        .saturating_mul(size_of::<ReferenceMinimizer>())
        .saturating_add(seed_hits.saturating_mul(size_of::<SeedHit>()))
        .saturating_add(unique_index_keys.saturating_mul(size_of::<MinimizerKey>()))
}

fn parse_max_memory_gb(value: &str) -> io::Result<u64> {
    let gb: f64 = value.parse::<f64>().map_err(|err| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("invalid --max-memory-gb value {value:?}: {err}"),
        )
    })?;
    if !gb.is_finite() || gb <= 0.0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--max-memory-gb must be a positive finite number",
        ));
    }

    Ok((gb * 1024.0 * 1024.0 * 1024.0) as u64)
}

#[cfg(debug_assertions)]
fn performance_metrics_enabled(args: &CliArgs) -> bool {
    args.verbose || env::var_os("FASTERANI_METRICS").is_some()
}

#[cfg(not(debug_assertions))]
fn performance_metrics_enabled(_args: &CliArgs) -> bool {
    false
}

struct QueryMappingStats {
    mapping_elapsed: std::time::Duration,
    summary_elapsed: std::time::Duration,
    mapping_count: usize,
    emitted_pairs: usize,
    #[cfg(debug_assertions)]
    max_mapping_result_bytes: usize,
    #[cfg(debug_assertions)]
    metrics: MappingMetrics,
}

struct RawQueryMappingStats {
    mapping_results: Vec<MappingResult>,
    mapping_elapsed: std::time::Duration,
    #[cfg(debug_assertions)]
    max_mapping_result_bytes: usize,
    #[cfg(debug_assertions)]
    metrics: MappingMetrics,
}

struct ShardQueryResult {
    shard_index: usize,
    reference_files: Vec<ReferenceFile>,
    reference_contig_names: Option<Vec<ReferenceContigName>>,
    mapping_results: Vec<MappingResult>,
    load_elapsed: std::time::Duration,
    mapping_elapsed: std::time::Duration,
    mapping_count: usize,
    #[cfg(debug_assertions)]
    max_mapping_result_bytes: usize,
    #[cfg(debug_assertions)]
    metrics: MappingMetrics,
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
    disable_reciprocal: bool,
) -> io::Result<(std::time::Duration, usize)> {
    let summary_start: Instant = Instant::now();
    let mapping_results_for_stats: Option<Vec<MappingResult>> = mapping_stats_output
        .as_ref()
        .map(|_| mapping_results.clone());
    let ani_computation: AniComputation = final_ani_computation(
        mapping_results,
        reference_files.len(),
        fragment_length,
        disable_reciprocal,
    );
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
        args.disable_reciprocal,
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
    let progress_enabled: bool = true;
    let performance_metrics_enabled: bool = performance_metrics_enabled(&args);
    let runtime_options: RuntimeOptions = RuntimeOptions {
        progress_enabled,
        max_memory_bytes: args.max_memory_bytes,
        worker_threads: args.threads,
    };
    check_memory_limit("startup", runtime_options)?;
    if progress_enabled
        && !fastani_compatible_fragment_mode(
            args.fragment_length,
            args.fragment_stride,
            args.min_fragment_length,
        )
    {
        eprintln!(
            "WARNING\tadaptive fragment mode enabled\tfragment_length={}\tfragment_stride={}\tmin_fragment_length={}",
            args.fragment_length, args.fragment_stride, args.min_fragment_length
        );
    }
    if progress_enabled {
        emit_progress(
            "parameters",
            &format!(
                "event=algorithm\tkmer_size={}\twindow_size={}\tfragment_length={}\tfragment_stride={}\tmin_fragment_length={}\tmin_identity={:.6}\tmash_confidence={:.6}\tminmer_count={}\tfreq_threshold_percent={:.6}\tsplit_n_run={}\tdisable_reciprocal={}",
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
                args.disable_reciprocal
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
        args.kmer_size,
        args.window_size,
        args.fragment_length,
        args.min_fragment_length,
        args.split_n_run,
        args.sketch_path.as_deref(),
        args.tmp_dir.as_deref(),
        args.bgzip,
        args.shard_size,
        args.shard_minimizers,
        args.index_build_mode,
        args.threads,
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
                                args.kmer_size,
                                args.window_size,
                                args.fragment_length,
                                args.min_fragment_length,
                                args.split_n_run,
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
                                shard_sketch.files.iter().cloned().collect();
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
                    args.disable_reciprocal,
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
