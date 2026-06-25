//! Seed-hit candidate discovery and sliding-window ANI scoring per query/reference pair.

use std::{cmp::Ordering, collections::HashSet, io};
#[cfg(debug_assertions)]
use std::{mem::size_of, time::Instant};

use noodles::fasta;
use rayon::prelude::*;

use crate::ani::{
    estimate_relaxed_minimum_shared_minimizers, query_fragment_ranges, query_fragment_sketch,
    split_sequence_ranges, AniComputation, AniSummary, MappingMetrics, MappingOutput,
    MappingResult, MappingResultKey, MappingScratch, QueryFile, QueryFragment, QueryFragmentSketch,
    ReferenceMinimizer, ReferenceSketch,
};
#[cfg(debug_assertions)]
use crate::ani::{MinimizerKey, QueryMemoryEstimate};

pub(crate) fn lower_bound_minimizer_position(
    minimizers: &[ReferenceMinimizer],
    position: u32,
) -> usize {
    minimizers.partition_point(|minimizer| minimizer.position < position)
}

impl QueryFile {
    /// Read a query FASTA file and split each contig into query fragments.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn collect(
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

    pub(crate) fn mapped_length(&self) -> u64 {
        self.mapped_length
    }

    pub(crate) fn total_minimizers(&self) -> usize {
        self.fragments
            .iter()
            .map(|fragment| fragment.minimizers.len())
            .sum()
    }

    pub(crate) fn total_seed_minimizers(&self) -> usize {
        self.fragments
            .iter()
            .map(|fragment| fragment.seed_minimizers.len())
            .sum()
    }

    #[cfg(debug_assertions)]
    pub(crate) fn memory_estimate(&self) -> QueryMemoryEstimate {
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
#[allow(clippy::too_many_arguments)]
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
#[allow(clippy::too_many_arguments)]
pub(crate) fn map_query_to_reference_parallel(
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
#[allow(clippy::too_many_arguments)]
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
        &mut scratch.slot_sorted_minimizers,
        &mut scratch.hit_ranges,
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
pub(crate) fn final_ani_computation(
    mut mapping_results: Vec<MappingResult>,
    reference_file_count: usize,
    fragment_length: u32,
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

    // Keep only reciprocal-best mappings: collapse query fragments that land in the
    // same reference position bin so each reference region is counted once.
    query_best_mappings.sort_by(|left, right| compare_refbin_bucket(left, right, fragment_length));

    let mut summary_mappings: Vec<MappingResult> = Vec::new();

    for mapping in query_best_mappings {
        if let Some(previous) = summary_mappings.last_mut() {
            if previous.reference_contig_id == mapping.reference_contig_id
                && reference_position_bin(previous.reference_start, fragment_length)
                    == reference_position_bin(mapping.reference_start, fragment_length)
            {
                *previous = mapping;
                continue;
            }
        }

        summary_mappings.push(mapping);
    }

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
