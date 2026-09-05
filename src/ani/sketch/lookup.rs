//! Query-time candidate-region discovery and scoring against a `ReferenceSketch`.

use crate::ani::{
    constants::MinimizerKey,
    mapping::lower_bound_minimizer_position,
    mash::{fastani_mash_distance, mash_distance_lower_bound},
    metrics::MappingMetrics,
    minimizer::SlidingSketchCounter,
    model::{
        query::{MappingResult, QueryFragment, ReferenceCandidateRegion},
        reference::{ReferenceMinimizer, ReferenceSketch, SeedHit},
    },
};

impl ReferenceSketch {
    fn database_frequency(&self, minimizer: &MinimizerKey, local_count: usize) -> usize {
        let Some(global_frequencies) = self.global_frequencies.as_ref() else {
            return local_count;
        };
        match global_frequencies.get(*minimizer) {
            Some(global_count) => global_count,
            None => {
                debug_assert!(
                    false,
                    "global-frequency artifact is missing a key present in a sketch shard"
                );
                local_count
            }
        }
    }

    /// Collect and merge seed-hit candidate intervals for one query fragment.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn find_candidate_regions(
        &self,
        query_minimizers: &[MinimizerKey],
        fragment_length: u32,
        minimum_shared_minimizers: usize,
        frequency_threshold: usize,
        seed_hits: &mut Vec<SeedHit>,
        candidate_regions: &mut Vec<ReferenceCandidateRegion>,
        slot_sorted_minimizers: &mut Vec<(u64, MinimizerKey)>,
        hit_ranges: &mut Vec<(usize, usize)>,
        #[cfg(debug_assertions)] mut mapping_metrics: Option<&mut MappingMetrics>,
    ) {
        seed_hits.clear();
        candidate_regions.clear();

        // Batch-compute MPHF slots and sort ascending so slot_keys / hit_offsets /
        // hit_counts are accessed sequentially rather than randomly across the mmap.
        self.index
            .slot_sorted_minimizers(query_minimizers, slot_sorted_minimizers);

        if slot_sorted_minimizers.is_empty() && !self.index.has_slot_sorted_lookup() {
            // Hash index path has no slot ordering; fall back to direct lookup.
            for minimizer in query_minimizers {
                let hits: Option<&[SeedHit]> = self.index.get(minimizer);
                let database_frequency: Option<usize> = hits
                    .map(<[SeedHit]>::len)
                    .map(|count| self.database_frequency(minimizer, count));
                #[cfg(debug_assertions)]
                if let Some(metrics) = mapping_metrics.as_deref_mut() {
                    metrics.record_seed_lookup(database_frequency, frequency_threshold);
                }
                if let Some(hits) = hits {
                    if database_frequency.expect("present hits have a frequency")
                        < frequency_threshold
                    {
                        seed_hits.extend_from_slice(hits);
                    }
                }
            }
        } else {
            hit_ranges.clear();
            hit_ranges.reserve(slot_sorted_minimizers.len());
            let mut accepted_hit_count: usize = 0;
            for &(slot, minimizer) in slot_sorted_minimizers.iter() {
                let hit_range = self.index.hit_range_by_slot(slot as usize, &minimizer);
                let database_frequency: Option<usize> = hit_range
                    .map(|(_, count)| count)
                    .map(|count| self.database_frequency(&minimizer, count));
                #[cfg(debug_assertions)]
                if let Some(metrics) = mapping_metrics.as_deref_mut() {
                    metrics.record_seed_lookup(database_frequency, frequency_threshold);
                }
                if let Some((offset, count)) = hit_range {
                    if database_frequency.expect("present hit range has a frequency")
                        < frequency_threshold
                    {
                        hit_ranges.push((offset, count));
                        accepted_hit_count = accepted_hit_count.saturating_add(count);
                    }
                }
            }

            hit_ranges.sort_unstable_by_key(|&(offset, _)| offset);
            seed_hits.reserve(accepted_hit_count);
            for &(offset, count) in hit_ranges.iter() {
                if let Some(hits) = self.index.hit_payload_range(offset, count) {
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

            let first_pos = first.position;
            let last_pos = last.position;

            if last_pos.saturating_sub(first_pos) >= fragment_length {
                continue;
            }

            let candidate_region = ReferenceCandidateRegion {
                reference_contig_id: first.reference_contig_id as usize,
                start_position: last_pos.saturating_sub(fragment_length - 1),
                end_position: first_pos,
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
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn score_candidate_region(
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
            if let Some(metrics) = mapping_metrics {
                metrics.scoring_window_steps += scoring_window_steps;
            }
        }

        if best_shared == 0 {
            return None;
        }

        let sketch_size: usize = query_fragment.minimizers.len();
        let jaccard: f64 = best_shared as f64 / sketch_size as f64;
        let distance: f64 = fastani_mash_distance(jaccard, kmer_size);
        let identity: f64 = 100.0 * (1.0 - distance);

        if min_identity > 0.0 {
            let lower_distance: f64 =
                mash_distance_lower_bound(distance, sketch_size, kmer_size, mash_confidence);
            let upper_identity: f64 = 100.0 * (1.0 - lower_distance);

            if upper_identity < min_identity {
                return None;
            }
        }

        Some(MappingResult {
            reference_file_id,
            reference_contig_id: candidate_region.reference_contig_id,
            query_fragment_id: query_fragment.id,
            query_fragment_length: query_fragment.length,
            reference_start: u32::midpoint(best_start, last_best_start),
            identity,
            query_minimizer_count: sketch_size,
            reference_minimizer_count: best_reference_minimizer_count,
            shared_minimizers: best_shared,
            union_minimizers: sketch_size,
            jaccard,
        })
    }
}
