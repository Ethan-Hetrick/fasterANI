//! Query-time candidate-region discovery and scoring against a `ReferenceSketch`.

use crate::ani::{
    fastani_mash_distance, lower_bound_minimizer_position, mash_distance_lower_bound,
    MappingMetrics, MappingResult, MinimizerKey, QueryFragment, ReferenceCandidateRegion,
    ReferenceMinimizer, ReferenceSketch, SeedHit, SlidingSketchCounter,
};

impl ReferenceSketch {
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

        seed_hits.sort_unstable_by_key(|hit| (hit.reference_contig_id, hit.minimizer_offset));

        // Helper: recover the genomic position for a hit from the contig minimizer array.
        let pos_of = |hit: SeedHit| -> u32 {
            self.contigs
                .minimizers(hit.reference_contig_id as usize)
                .and_then(|mins| mins.get(hit.minimizer_offset as usize))
                .map_or(0, |m| m.position)
        };

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

            let first_pos = pos_of(first);
            let last_pos  = pos_of(last);

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
