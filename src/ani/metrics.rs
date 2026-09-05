//! Hot-path mapping counters and optional debug-only seed-hit histograms.

use crate::ani::model::query::MappingResult;

/// Fine-grained debug-only seed-hit histogram buckets.
#[cfg(debug_assertions)]
pub(crate) const SEED_HIT_HISTOGRAM_UPPER_BOUNDS: [usize; 13] =
    [1, 2, 4, 9, 24, 49, 99, 249, 499, 999, 4_999, 9_999, 49_999];

#[cfg(debug_assertions)]
pub(crate) const SEED_HIT_HISTOGRAM_OVERFLOW_LABEL: &str = "50000+";

#[derive(Clone, Default)]
pub(crate) struct MappingMetrics {
    pub(crate) candidate_regions_found: usize,
    pub(crate) candidate_regions_scored: usize,
    #[cfg(debug_assertions)]
    pub(crate) candidate_discovery_calls: usize,
    #[cfg(debug_assertions)]
    pub(crate) seed_hits_collected: usize,
    #[cfg(debug_assertions)]
    pub(crate) reference_minimizers_scanned: usize,
    #[cfg(debug_assertions)]
    pub(crate) scoring_window_steps: usize,
    #[cfg(debug_assertions)]
    pub(crate) retained_mappings: usize,
    pub(crate) candidate_discovery_elapsed: std::time::Duration,
    #[cfg(debug_assertions)]
    pub(crate) scoring_elapsed: std::time::Duration,
    #[cfg(debug_assertions)]
    pub(crate) seed_lookup_count: usize,
    #[cfg(debug_assertions)]
    pub(crate) seed_lookup_zero_hits: usize,
    #[cfg(debug_assertions)]
    pub(crate) seed_lookup_skipped_by_frequency: usize,
    #[cfg(debug_assertions)]
    pub(crate) seed_hit_list_max: usize,
    #[cfg(debug_assertions)]
    pub(crate) seed_hit_list_bins: [usize; SEED_HIT_HISTOGRAM_UPPER_BOUNDS.len() + 1],
    #[cfg(debug_assertions)]
    pub(crate) seed_hit_list_bin_hits: [usize; SEED_HIT_HISTOGRAM_UPPER_BOUNDS.len() + 1],
}

impl MappingMetrics {
    pub(crate) fn merge(&mut self, other: Self) {
        self.candidate_regions_found += other.candidate_regions_found;
        self.candidate_regions_scored += other.candidate_regions_scored;
        self.candidate_discovery_elapsed += other.candidate_discovery_elapsed;
        #[cfg(debug_assertions)]
        {
            self.candidate_discovery_calls += other.candidate_discovery_calls;
            self.seed_hits_collected += other.seed_hits_collected;
            self.reference_minimizers_scanned += other.reference_minimizers_scanned;
            self.scoring_window_steps += other.scoring_window_steps;
            self.retained_mappings += other.retained_mappings;
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
    }

    #[cfg(debug_assertions)]
    pub(crate) fn record_seed_lookup(
        &mut self,
        hit_list_len: Option<usize>,
        frequency_threshold: usize,
    ) {
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

#[cfg(debug_assertions)]
fn seed_hit_histogram_bin_index(hit_list_len: usize) -> usize {
    SEED_HIT_HISTOGRAM_UPPER_BOUNDS
        .iter()
        .position(|&upper_bound| hit_list_len <= upper_bound)
        .unwrap_or(SEED_HIT_HISTOGRAM_UPPER_BOUNDS.len())
}

/// Raw mapping results plus optional hot-path metrics for one query file.
pub(crate) struct MappingOutput {
    pub(crate) results: Vec<MappingResult>,
    pub(crate) metrics: MappingMetrics,
}
