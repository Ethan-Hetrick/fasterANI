//! Reference-side data types: files, contigs, minimizers, indexes, and sketch metadata.

use std::cmp::Reverse;
use std::collections::HashMap;
use std::path::Path;

use boomphf::Mphf;
use serde::{Deserialize, Serialize};

use crate::ani::{
    default_fragment_length, default_shard_minimizers, IndexBuildMode, MinimizerKey,
    MmapReferenceContigs, MmapReferenceIndex, ReferenceHitMap,
};

/// The FastANI-style algorithm parameters that travel together through sketch
/// construction, persistence, and loading. Grouping them keeps the long build
/// and load signatures readable and impossible to mis-order.
#[derive(Clone, Copy, Debug)]
pub(crate) struct SketchParams {
    pub(crate) kmer_size: usize,
    pub(crate) window_size: usize,
    pub(crate) fragment_length: u32,
    pub(crate) min_fragment_length: u32,
    pub(crate) split_n_run: usize,
}

/// Filesystem and sharding options used when building a sharded reference database.
#[derive(Clone, Copy)]
pub(crate) struct ShardedBuildOptions<'a> {
    pub(crate) tmp_dir: Option<&'a Path>,
    pub(crate) bgzip: bool,
    pub(crate) shard_size: usize,
    pub(crate) shard_minimizers: usize,
    pub(crate) index_build_mode: IndexBuildMode,
    pub(crate) threads: usize,
}

/// One input reference genome file and the number of bases FastANI considers mappable.
#[derive(Clone, Deserialize, Serialize)]
pub(crate) struct ReferenceFile {
    pub(crate) path: String,
    pub(crate) mapped_length: u64,
}

/// One FASTA record from a reference file and its ordered minimizers.
#[derive(Clone, Deserialize, Serialize)]
pub(crate) struct ReferenceContig {
    pub(crate) file_id: usize,
    pub(crate) minimizers: Vec<ReferenceMinimizer>,
}

/// Optional human-readable metadata for one indexed reference contig or split segment.
#[derive(Clone)]
pub(crate) struct ReferenceContigName {
    pub(crate) file_id: usize,
    pub(crate) name: String,
    pub(crate) segment_start: u32,
    pub(crate) segment_end: u32,
}

/// Compact persisted metadata for one reference contig in a memory-mapped sketch.
#[repr(C)]
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct ContigRecord {
    pub(crate) minimizer_offset: u64,
    pub(crate) file_id: u32,
    pub(crate) minimizer_count: u32,
}

/// A reference minimizer hash and its zero-based contig position.
#[repr(C)]
#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub(crate) struct ReferenceMinimizer {
    pub(crate) hash: MinimizerKey,
    pub(crate) position: u32,
}

/// Compact seed hit stored in the reference index.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, Deserialize, Eq, PartialEq, Serialize)]
pub(crate) struct SeedHit {
    pub(crate) reference_contig_id: u32,
    pub(crate) position: u32,
}

/// Fully loaded reference sketch, backed by either an in-memory hash map or an mmap cache.
pub(crate) struct ReferenceSketch {
    pub(crate) files: Vec<ReferenceFile>,
    pub(crate) contigs: ReferenceContigs,
    pub(crate) contig_names: Option<Vec<ReferenceContigName>>,
    pub(crate) index: ReferenceIndex,
}

/// Reference contig minimizers, either owned after a fresh build or sliced from a loaded sketch.
pub(crate) enum ReferenceContigs {
    Owned(Vec<ReferenceContig>),
    Mmap(MmapReferenceContigs),
}

/// Lookup table from minimizer hash to all reference positions containing that minimizer.
pub(crate) enum ReferenceIndex {
    Hash(ReferenceHitMap),
    Mphf(MmapReferenceIndex),
}

/// JSON metadata stored at the front of a `.fasketch` cache.
#[derive(Deserialize, Serialize)]
pub(crate) struct CachedReferenceMetadata {
    pub(crate) version: u32,
    pub(crate) k: usize,
    pub(crate) w: usize,
    #[serde(default)]
    pub(crate) key_mode: String,
    #[serde(default = "default_fragment_length")]
    pub(crate) fragment_length: u32,
    pub(crate) min_fragment_length: u32,
    #[serde(default)]
    pub(crate) split_n_run: usize,
    #[serde(default)]
    pub(crate) dust_enabled: bool,
    pub(crate) files: Vec<ReferenceFile>,
    pub(crate) mphf: Mphf<MinimizerKey>,
    pub(crate) key_count: usize,
    pub(crate) hit_count: usize,
    pub(crate) contig_count: usize,
    pub(crate) reference_minimizer_count: usize,
}

/// Top-level metadata for a manifest-backed sharded reference database.
#[derive(Clone, Deserialize, Serialize)]
pub(crate) struct ShardManifest {
    pub(crate) sketch_format_version: u32,
    pub(crate) database_schema_version: u32,
    pub(crate) k: usize,
    pub(crate) w: usize,
    pub(crate) key_mode: String,
    #[serde(default = "default_fragment_length")]
    pub(crate) fragment_length: u32,
    pub(crate) min_fragment_length: u32,
    pub(crate) split_n_run: usize,
    #[serde(default)]
    pub(crate) dust_enabled: bool,
    pub(crate) shard_size: usize,
    #[serde(default = "default_shard_minimizers")]
    pub(crate) shard_minimizers: usize,
    pub(crate) total_references: usize,
    pub(crate) total_reference_contigs: usize,
    pub(crate) total_mapped_reference_length: u64,
    pub(crate) total_reference_minimizers: usize,
    pub(crate) total_shard_unique_minimizers: usize,
    pub(crate) build_unix_seconds: u64,
    pub(crate) build_args: Vec<String>,
    pub(crate) reference_list_checksum: u64,
    pub(crate) shards: Vec<ShardManifestEntry>,
}

/// One physical sketch shard referenced by a sharded database manifest.
#[derive(Clone, Deserialize, Serialize)]
pub(crate) struct ShardManifestEntry {
    pub(crate) shard_index: usize,
    pub(crate) filename: String,
    pub(crate) first_reference: usize,
    pub(crate) reference_count: usize,
    pub(crate) reference_contigs: usize,
    pub(crate) mapped_reference_length: u64,
    pub(crate) reference_minimizers: usize,
    pub(crate) unique_minimizers: usize,
}

/// Small summary returned by a streaming shard build.
#[derive(Clone, Copy, Default)]
pub(crate) struct SketchBuildStats {
    pub(crate) reference_count: usize,
    pub(crate) reference_contig_count: usize,
    pub(crate) mapped_reference_length: u64,
    pub(crate) reference_minimizer_count: usize,
    pub(crate) unique_minimizer_count: usize,
}

/// Completed shard build result used to assemble the final database manifest.
pub(crate) struct ShardBuildResult {
    pub(crate) entry: ShardManifestEntry,
}

/// Reference-list range selected for one physical sketch shard.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(crate) struct ShardPlan {
    pub(crate) first_reference: usize,
    pub(crate) reference_count: usize,
    pub(crate) estimated_minimizers: usize,
}

#[cfg(debug_assertions)]
#[allow(dead_code)]
#[derive(Default)]
pub(crate) struct ReferenceMemoryEstimate {
    pub(crate) reference_minimizers: usize,
    pub(crate) reference_minimizer_vec_bytes: usize,
    pub(crate) unique_index_keys: usize,
    pub(crate) seed_hits: usize,
    pub(crate) seed_hit_vec_bytes: usize,
    pub(crate) hash_index_rough_bytes: usize,
    pub(crate) mmap_file_bytes: usize,
    pub(crate) mmap_slot_key_bytes: usize,
    pub(crate) mmap_hit_offset_bytes: usize,
    pub(crate) mmap_hit_count_bytes: usize,
    pub(crate) mmap_hit_payload_bytes: usize,
    pub(crate) mmap_contig_record_bytes: usize,
    pub(crate) mmap_reference_minimizer_bytes: usize,
}

impl ReferenceIndex {
    pub(crate) fn get(&self, minimizer: &MinimizerKey) -> Option<&[SeedHit]> {
        match self {
            Self::Hash(index) => index.get(minimizer).map(Vec::as_slice),
            Self::Mphf(index) => index.get(minimizer),
        }
    }

    /// Populate `out` with `(slot, minimizer)` pairs for all minimizers that
    /// exist in the index, sorted ascending by slot. Clears `out` first.
    pub(crate) fn slot_sorted_minimizers(
        &self,
        minimizers: &[MinimizerKey],
        out: &mut Vec<(u64, MinimizerKey)>,
    ) {
        out.clear();
        match self {
            Self::Mphf(index) => {
                for &minimizer in minimizers {
                    if let Some(slot) = index.mphf.try_hash(&minimizer) {
                        if (slot as usize) < index.key_count {
                            out.push((slot, minimizer));
                        }
                    }
                }
                out.sort_unstable_by_key(|&(slot, _)| slot);
            }
            Self::Hash(_) => {
                // In-memory hash map has no meaningful slot ordering;
                // leave out empty and fall back to direct lookup in the caller.
            }
        }
    }

    pub(crate) fn get_by_slot(&self, slot: usize, minimizer: &MinimizerKey) -> Option<&[SeedHit]> {
        match self {
            Self::Mphf(index) => index.get_by_slot(slot, minimizer),
            Self::Hash(_) => self.get(minimizer),
        }
    }

    pub(crate) fn len(&self) -> usize {
        match self {
            Self::Hash(index) => index.len(),
            Self::Mphf(index) => index.key_count,
        }
    }

    pub(crate) fn frequency_threshold(&self, percent: f64) -> usize {
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
        frequencies.sort_unstable_by_key(|&(frequency, _count)| Reverse(frequency));

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
