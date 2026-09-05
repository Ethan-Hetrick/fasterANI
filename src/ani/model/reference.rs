//! Reference-side data types: files, contigs, minimizers, indexes, and sketch metadata.

use std::cmp::{Ordering, Reverse};
use std::collections::HashMap;
use std::{path::Path, sync::Arc};

use boomphf::Mphf;
use serde::{Deserialize, Serialize};

use crate::ani::{
    constants::{MinimizerKey, ReferenceHitMap},
    mmap::{MmapReferenceContigs, MmapReferenceIndex},
    sketch::frequency::GlobalFrequencyIndex,
    validation::{default_fragment_length, default_max_shard_size_bytes},
};

/// The FastANI-style algorithm parameters that travel together through sketch
/// construction, persistence, and loading. Grouping them keeps the long build
/// and load signatures readable and impossible to mis-order.
#[derive(Clone, Copy, Debug)]
pub(crate) struct SketchParams {
    pub(crate) kmer_size: usize,
    pub(crate) window_size: usize,
    pub(crate) minimizer_hash_seed: u32,
    pub(crate) fragment_length: u32,
    pub(crate) min_fragment_length: u32,
    pub(crate) split_n_run: usize,
}

/// Filesystem and sharding options used when building a sharded reference database.
#[derive(Clone, Copy)]
pub(crate) struct ShardedBuildOptions<'a> {
    pub(crate) tmp_dir: Option<&'a Path>,
    pub(crate) max_shard_size_bytes: u64,
    pub(crate) threads: usize,
    pub(crate) force_rebuild: bool,
}

/// One input reference genome file and the number of bases `FastANI` considers mappable.
#[derive(Clone, Deserialize, Serialize)]
pub(crate) struct ReferenceFile {
    pub(crate) path: String,
    pub(crate) mapped_length: u64,
    #[serde(default)]
    pub(crate) original_length: u64,
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

const _: [(); 16] = [(); std::mem::size_of::<ContigRecord>()];
const _: [(); 0] = [(); std::mem::offset_of!(ContigRecord, minimizer_offset)];
const _: [(); 8] = [(); std::mem::offset_of!(ContigRecord, file_id)];
const _: [(); 12] = [(); std::mem::offset_of!(ContigRecord, minimizer_count)];

/// A reference minimizer hash and its zero-based contig position.
#[repr(C)]
#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub(crate) struct ReferenceMinimizer {
    pub(crate) hash: MinimizerKey,
    pub(crate) position: u32,
}

const _: [(); 8] = [(); std::mem::size_of::<ReferenceMinimizer>()];
const _: [(); 0] = [(); std::mem::offset_of!(ReferenceMinimizer, hash)];
const _: [(); 4] = [(); std::mem::offset_of!(ReferenceMinimizer, position)];

/// Compact seed hit stored in the reference index.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, Deserialize, Eq, PartialEq, Serialize)]
pub(crate) struct SeedHit {
    pub(crate) reference_contig_id: u32,
    pub(crate) position: u32,
}

const _: [(); 8] = [(); std::mem::size_of::<SeedHit>()];
const _: [(); 0] = [(); std::mem::offset_of!(SeedHit, reference_contig_id)];
const _: [(); 4] = [(); std::mem::offset_of!(SeedHit, position)];

/// Fully loaded reference sketch, backed by either an in-memory hash map or an mmap cache.
pub(crate) struct ReferenceSketch {
    pub(crate) files: Vec<ReferenceFile>,
    pub(crate) contigs: ReferenceContigs,
    pub(crate) contig_names: Option<Vec<ReferenceContigName>>,
    pub(crate) index: ReferenceIndex,
    pub(crate) global_frequencies: Option<Arc<GlobalFrequencyIndex>>,
}

/// Reference contig minimizers, either owned after a fresh build or sliced from a loaded sketch.
pub(crate) enum ReferenceContigs {
    Owned(Vec<ReferenceContig>),
    Mmap(MmapReferenceContigs),
}

/// Lookup table from minimizer hash to all reference positions containing that minimizer.
#[allow(dead_code)]
pub(crate) enum ReferenceIndex {
    Transient(TransientReferenceIndex),
    Hash(ReferenceHitMap),
    HashShards(Vec<ReferenceHashShard>),
    Mphf(MmapReferenceIndex),
}

/// Flat one-shot lookup table for freshly built no-save reference sketches.
pub(crate) struct TransientReferenceIndex {
    pub(crate) keys: Vec<MinimizerKey>,
    pub(crate) hit_offsets: Vec<usize>,
    pub(crate) hit_payloads: Vec<SeedHit>,
}

/// One key-range shard of an in-memory reference hit map.
pub(crate) struct ReferenceHashShard {
    pub(crate) first_key: MinimizerKey,
    pub(crate) last_key: MinimizerKey,
    pub(crate) index: ReferenceHitMap,
}

/// JSON metadata stored at the front of a `.fasketch` cache.
#[derive(Deserialize, Serialize)]
pub(crate) struct CachedReferenceMetadata {
    pub(crate) version: u32,
    pub(crate) k: usize,
    pub(crate) w: usize,
    #[serde(default)]
    pub(crate) minimizer_hash_seed: u32,
    #[serde(default = "default_fragment_length")]
    pub(crate) fragment_length: u32,
    pub(crate) min_fragment_length: u32,
    #[serde(default)]
    pub(crate) split_n_run: usize,
    pub(crate) reference_count: usize,
    pub(crate) mphf: Mphf<MinimizerKey>,
    pub(crate) key_count: usize,
    pub(crate) hit_count: usize,
    pub(crate) contig_count: usize,
    pub(crate) reference_minimizer_count: usize,
    pub(crate) name_sidecar_filename: String,
    pub(crate) name_sidecar_file_bytes: u64,
}

/// Top-level metadata for a manifest-backed sharded reference database.
#[derive(Clone, Deserialize, Serialize)]
pub(crate) struct ShardManifest {
    pub(crate) sketch_format_version: u32,
    pub(crate) database_schema_version: u32,
    pub(crate) k: usize,
    pub(crate) w: usize,
    #[serde(default)]
    pub(crate) minimizer_hash_seed: u32,
    #[serde(default = "default_fragment_length")]
    pub(crate) fragment_length: u32,
    pub(crate) min_fragment_length: u32,
    pub(crate) split_n_run: usize,
    #[serde(default = "default_max_shard_size_bytes")]
    pub(crate) max_shard_size_bytes: u64,
    pub(crate) total_references: usize,
    pub(crate) total_reference_contigs: usize,
    pub(crate) total_mapped_reference_length: u64,
    pub(crate) total_reference_minimizers: usize,
    pub(crate) total_shard_unique_minimizers: usize,
    pub(crate) build_unix_seconds: u64,
    #[serde(default)]
    pub(crate) generation_id: String,
    #[serde(default)]
    pub(crate) global_frequency_filename: String,
    #[serde(default)]
    pub(crate) global_frequency_file_bytes: u64,
    #[serde(default)]
    pub(crate) total_unique_minimizers: usize,
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
    #[serde(default)]
    pub(crate) estimated_file_bytes: u64,
    #[serde(default)]
    pub(crate) file_bytes: u64,
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
    pub(crate) estimated_file_bytes: u64,
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

impl TransientReferenceIndex {
    pub(crate) fn from_sorted_hits_with_key_count(
        hits: &[(MinimizerKey, SeedHit)],
        unique_keys: usize,
    ) -> Self {
        let mut keys: Vec<MinimizerKey> = Vec::with_capacity(unique_keys);
        let mut hit_offsets: Vec<usize> = Vec::with_capacity(unique_keys.saturating_add(1));
        let mut hit_payloads: Vec<SeedHit> = Vec::with_capacity(hits.len());
        let mut cursor: usize = 0usize;

        while cursor < hits.len() {
            let key: MinimizerKey = hits[cursor].0;
            keys.push(key);
            hit_offsets.push(hit_payloads.len());

            while cursor < hits.len() && hits[cursor].0 == key {
                hit_payloads.push(hits[cursor].1);
                cursor += 1;
            }
        }
        hit_offsets.push(hit_payloads.len());

        Self {
            keys,
            hit_offsets,
            hit_payloads,
        }
    }

    pub(crate) fn get(&self, minimizer: &MinimizerKey) -> Option<&[SeedHit]> {
        let slot: usize = self.keys.binary_search(minimizer).ok()?;
        let start: usize = self.hit_offsets[slot];
        let end: usize = self.hit_offsets[slot + 1];
        self.hit_payloads.get(start..end)
    }

    pub(crate) fn hit_range_by_slot(
        &self,
        slot: usize,
        minimizer: &MinimizerKey,
    ) -> Option<(usize, usize)> {
        if self.keys.get(slot)? != minimizer {
            return None;
        }
        let start: usize = self.hit_offsets[slot];
        let end: usize = self.hit_offsets[slot + 1];
        Some((start, end - start))
    }

    pub(crate) fn hit_payload_range(&self, offset: usize, count: usize) -> Option<&[SeedHit]> {
        let end: usize = offset.checked_add(count)?;
        self.hit_payloads.get(offset..end)
    }
}

impl ReferenceIndex {
    pub(crate) fn get(&self, minimizer: &MinimizerKey) -> Option<&[SeedHit]> {
        match self {
            Self::Transient(index) => index.get(minimizer),
            Self::Hash(index) => index.get(minimizer).map(Vec::as_slice),
            Self::HashShards(shards) => {
                let shard_index: usize =
                    shards.partition_point(|shard| shard.last_key < *minimizer);
                let shard: &ReferenceHashShard = shards.get(shard_index)?;
                if *minimizer < shard.first_key {
                    return None;
                }
                shard.index.get(minimizer).map(Vec::as_slice)
            }
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
            Self::Transient(index) => {
                for &minimizer in minimizers {
                    if let Ok(slot) = index.keys.binary_search(&minimizer) {
                        out.push((slot as u64, minimizer));
                    }
                }
                out.sort_unstable_by_key(|&(slot, _)| slot);
            }
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
            Self::Hash(_) | Self::HashShards(_) => {
                // In-memory hash map has no meaningful slot ordering;
                // leave out empty and fall back to direct lookup in the caller.
            }
        }
    }

    pub(crate) fn has_slot_sorted_lookup(&self) -> bool {
        matches!(self, Self::Transient(_) | Self::Mphf(_))
    }

    pub(crate) fn hit_range_by_slot(
        &self,
        slot: usize,
        minimizer: &MinimizerKey,
    ) -> Option<(usize, usize)> {
        match self {
            Self::Transient(index) => index.hit_range_by_slot(slot, minimizer),
            Self::Mphf(index) => index
                .hit_range_by_slot(slot, minimizer)
                .map(|(offset, count)| (offset as usize, count as usize)),
            Self::Hash(_) | Self::HashShards(_) => None,
        }
    }

    pub(crate) fn hit_payload_range(&self, offset: usize, count: usize) -> Option<&[SeedHit]> {
        match self {
            Self::Transient(index) => index.hit_payload_range(offset, count),
            Self::Mphf(index) => {
                let offset = u32::try_from(offset).ok()?;
                let count = u32::try_from(count).ok()?;
                index.hit_payload_range(offset, count)
            }
            Self::Hash(_) | Self::HashShards(_) => None,
        }
    }

    pub(crate) fn len(&self) -> usize {
        match self {
            Self::Transient(index) => index.keys.len(),
            Self::Hash(index) => index.len(),
            Self::HashShards(shards) => shards.iter().map(|shard| shard.index.len()).sum(),
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
            Self::Transient(index) => {
                for hit_count in index
                    .hit_offsets
                    .windows(2)
                    .map(|offsets| offsets[1] - offsets[0])
                {
                    *histogram.entry(hit_count).or_default() += 1;
                    total_unique_minimizers += 1;
                }
            }
            Self::Hash(index) => {
                for hits in index.values() {
                    *histogram.entry(hits.len()).or_default() += 1;
                    total_unique_minimizers += 1;
                }
            }
            Self::HashShards(shards) => {
                for shard in shards {
                    for hits in shard.index.values() {
                        *histogram.entry(hits.len()).or_default() += 1;
                        total_unique_minimizers += 1;
                    }
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
            match sum.cmp(&minimizers_to_ignore) {
                Ordering::Less => threshold = frequency,
                Ordering::Equal => {
                    threshold = frequency;
                    break;
                }
                Ordering::Greater => break,
            }
        }

        threshold
    }
}

#[cfg(test)]
mod tests {
    use super::{MinimizerKey, ReferenceHitMap, ReferenceIndex, SeedHit, TransientReferenceIndex};

    fn hit(reference_contig_id: u32, position: u32) -> SeedHit {
        SeedHit {
            reference_contig_id,
            position,
        }
    }

    #[test]
    fn transient_reference_index_returns_present_key_hits() {
        let hits: Vec<(MinimizerKey, SeedHit)> = vec![
            (10, hit(0, 1)),
            (10, hit(0, 3)),
            (25, hit(1, 8)),
            (40, hit(2, 13)),
        ];
        let index = TransientReferenceIndex::from_sorted_hits_with_key_count(&hits, 3);

        assert_eq!(index.get(&10), Some([hit(0, 1), hit(0, 3)].as_slice()));
        assert_eq!(index.get(&25), Some([hit(1, 8)].as_slice()));
    }

    #[test]
    fn transient_reference_index_returns_none_for_missing_key() {
        let hits: Vec<(MinimizerKey, SeedHit)> = vec![(10, hit(0, 1)), (25, hit(1, 8))];
        let index = TransientReferenceIndex::from_sorted_hits_with_key_count(&hits, 2);

        assert_eq!(index.get(&11), None);
    }

    #[test]
    fn transient_reference_index_groups_repeated_keys_into_offsets() {
        let hits: Vec<(MinimizerKey, SeedHit)> = vec![
            (7, hit(0, 1)),
            (7, hit(0, 2)),
            (7, hit(0, 3)),
            (9, hit(1, 5)),
        ];
        let index = TransientReferenceIndex::from_sorted_hits_with_key_count(&hits, 2);

        assert_eq!(index.keys, vec![7, 9]);
        assert_eq!(index.hit_offsets, vec![0, 3, 4]);
        assert_eq!(
            index.hit_payloads,
            vec![hit(0, 1), hit(0, 2), hit(0, 3), hit(1, 5)]
        );
    }

    #[test]
    fn transient_reference_index_frequency_threshold_matches_hash() {
        let hits: Vec<(MinimizerKey, SeedHit)> = vec![
            (5, hit(0, 1)),
            (5, hit(0, 2)),
            (5, hit(0, 3)),
            (8, hit(1, 1)),
            (8, hit(1, 2)),
            (13, hit(2, 1)),
        ];
        let transient = ReferenceIndex::Transient(
            TransientReferenceIndex::from_sorted_hits_with_key_count(&hits, 3),
        );
        let mut hash: ReferenceHitMap = ReferenceHitMap::default();
        for (key, hit) in hits {
            hash.entry(key).or_default().push(hit);
        }
        let hash = ReferenceIndex::Hash(hash);

        assert_eq!(
            transient.frequency_threshold(50.0),
            hash.frequency_threshold(50.0)
        );
        assert_eq!(
            transient.frequency_threshold(100.0),
            hash.frequency_threshold(100.0)
        );
    }
}
