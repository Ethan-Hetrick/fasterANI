//! Query-side data types: fragments, mapping results, and ANI summaries.

use std::collections::HashSet;

use crate::ani::{MinimizerKey, SeedHit, SlidingSketchCounter};

/// One fixed-length query fragment with full scoring minimizers and optional reduced seed hashes.
#[derive(Clone)]
pub(crate) struct QueryFragment {
    pub(crate) id: usize,
    pub(crate) contig_id: usize,
    pub(crate) start: u32,
    pub(crate) end: u32,
    pub(crate) length: u32,
    pub(crate) minimizers: Vec<MinimizerKey>,
    pub(crate) seed_minimizers: Vec<MinimizerKey>,
}

/// Query genome split into FastANI-style fixed-width fragments.
pub(crate) struct QueryFile {
    pub(crate) fragments: Vec<QueryFragment>,
    pub(crate) contig_names: Vec<String>,
    pub(crate) mapped_length: u64,
}

/// Candidate reference position interval produced by clustered seed hits.
#[derive(Clone, Copy)]
pub(crate) struct ReferenceCandidateRegion {
    pub(crate) reference_contig_id: usize,
    pub(crate) start_position: u32,
    pub(crate) end_position: u32,
}

/// One retained mapping between a query fragment and a reference location.
#[derive(Clone)]
pub(crate) struct MappingResult {
    pub(crate) reference_file_id: usize,
    pub(crate) reference_contig_id: usize,
    pub(crate) query_fragment_id: usize,
    pub(crate) query_fragment_length: u32,
    pub(crate) reference_start: u32,
    pub(crate) identity: f64,
    pub(crate) query_minimizer_count: usize,
    pub(crate) reference_minimizer_count: usize,
    pub(crate) shared_minimizers: usize,
    pub(crate) union_minimizers: usize,
    pub(crate) jaccard: f64,
}

/// Stable identity for marking mappings that survive reciprocal-best filtering.
#[derive(Hash, Eq, PartialEq)]
pub(crate) struct MappingResultKey {
    pub(crate) reference_file_id: usize,
    pub(crate) reference_contig_id: usize,
    pub(crate) query_fragment_id: usize,
    pub(crate) query_fragment_length: u32,
    pub(crate) reference_start: u32,
    pub(crate) identity_bits: u64,
    pub(crate) query_minimizer_count: usize,
    pub(crate) reference_minimizer_count: usize,
    pub(crate) shared_minimizers: usize,
    pub(crate) union_minimizers: usize,
}

impl MappingResultKey {
    pub(crate) fn from_mapping(mapping: &MappingResult) -> Self {
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
pub(crate) struct AniSummary {
    pub(crate) shared_fragments: usize,
    pub(crate) shared_bases: u64,
    pub(crate) weighted_identity_sum: f64,
    pub(crate) distribution_stats: AniDistributionStats,
}

/// Distribution statistics for retained fragment ANI values.
#[derive(Clone, Copy)]
pub struct AniDistributionStats {
    pub median: f64,
    pub stddev: f64,
    pub ci_95_lower: f64,
    pub ci_95_upper: f64,
    pub p99: f64,
    pub p80: f64,
}

impl Default for AniDistributionStats {
    fn default() -> Self {
        Self {
            median: f64::NAN,
            stddev: f64::NAN,
            ci_95_lower: f64::NAN,
            ci_95_upper: f64::NAN,
            p99: f64::NAN,
            p80: f64::NAN,
        }
    }
}

/// Final per-reference summaries plus the exact mappings that contributed to them.
pub(crate) struct AniComputation {
    pub(crate) summaries: Vec<AniSummary>,
    pub(crate) reciprocal_best_keys: HashSet<MappingResultKey>,
}

/// Per-thread reusable buffers for query-to-reference mapping.
#[derive(Default)]
pub(crate) struct MappingScratch {
    pub(crate) seed_hits: Vec<SeedHit>,
    pub(crate) candidate_regions: Vec<ReferenceCandidateRegion>,
    pub(crate) fragment_mappings: Vec<MappingResult>,
    pub(crate) counter: SlidingSketchCounter,
    pub(crate) slot_sorted_minimizers: Vec<(u64, MinimizerKey)>,
    pub(crate) hit_ranges: Vec<(usize, usize)>,
}
