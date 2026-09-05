//! Seed-hit candidate discovery and sliding-window ANI scoring per query/reference pair.

use std::{cmp::Ordering, collections::HashSet, io, mem::size_of, time::Instant};

use noodles::fasta;
use rayon::prelude::*;

use crate::ani::{
    constants::MinimizerKey,
    mash::{binomial_survival, estimate_relaxed_minimum_shared_minimizers, t_cdf_approx},
    metrics::{MappingMetrics, MappingOutput},
    minimizer::{
        query_fragment_ranges, query_fragment_sketch, split_sequence_ranges, QueryFragmentSketch,
    },
    model::{
        AniComputation, AniDistributionStats, AniSummary, ContigAniSummary, MappingResult,
        MappingResultKey, MappingScratch, QueryFile, QueryFragment, ReferenceMinimizer,
        ReferenceSketch,
    },
};

/// Run-scoped executor reused by every query/reference mapping operation.
pub(crate) struct MappingExecutor {
    pool: Option<rayon::ThreadPool>,
}

impl MappingExecutor {
    pub(crate) fn new(threads: usize) -> io::Result<Self> {
        if threads <= 1 {
            return Ok(Self { pool: None });
        }

        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .map_err(|err| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("failed to initialize mapping thread pool: {err}"),
                )
            })?;

        Ok(Self { pool: Some(pool) })
    }
}

fn checked_query_coordinate(value: usize, field: &str, contig_name: &str) -> io::Result<u32> {
    u32::try_from(value).map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "query contig '{contig_name}' {field} coordinate {value} exceeds the u32 limit"
            ),
        )
    })
}

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
        minimizer_hash_seed: u32,
        minmer_count: Option<usize>,
        fragment_length: u32,
        fragment_stride: u32,
        min_fragment_length: u32,
        split_n_run: usize,
        allow_empty_fragments: bool,
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
                        minimizer_hash_seed,
                        minmer_count,
                    );
                    if fragment_sketch.minimizers.is_empty() {
                        continue;
                    }

                    let contig_name: &str = &contig_names[contig_id];
                    let fragment_length: u32 = checked_query_coordinate(
                        fragment_range.len(),
                        "fragment length",
                        contig_name,
                    )?;
                    let query_start: u32 = checked_query_coordinate(
                        segment_start + fragment_range.start,
                        "fragment start",
                        contig_name,
                    )?;
                    let query_end: u32 = checked_query_coordinate(
                        segment_start + fragment_range.end,
                        "fragment end",
                        contig_name,
                    )?;
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

        if fragments.is_empty() && !allow_empty_fragments {
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

    pub(crate) fn estimated_owned_bytes(&self) -> usize {
        let fragment_struct_bytes: usize = self
            .fragments
            .capacity()
            .saturating_mul(size_of::<QueryFragment>());
        let query_minimizer_vec_bytes: usize =
            self.fragments.iter().fold(0usize, |total, fragment| {
                total.saturating_add(
                    fragment
                        .minimizers
                        .capacity()
                        .saturating_mul(size_of::<MinimizerKey>()),
                )
            });
        let seed_minimizer_vec_bytes: usize =
            self.fragments.iter().fold(0usize, |total, fragment| {
                total.saturating_add(
                    fragment
                        .seed_minimizers
                        .capacity()
                        .saturating_mul(size_of::<MinimizerKey>()),
                )
            });
        let contig_name_bytes: usize = self
            .contig_names
            .capacity()
            .saturating_mul(size_of::<String>())
            .saturating_add(
                self.contig_names
                    .iter()
                    .fold(0usize, |total, name| total.saturating_add(name.capacity())),
            );
        fragment_struct_bytes
            .saturating_add(query_minimizer_vec_bytes)
            .saturating_add(seed_minimizer_vec_bytes)
            .saturating_add(contig_name_bytes)
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
        metrics: mapping_metrics,
    }
}

/// Map one query using the Rayon pool that is already active on the current thread.
#[allow(clippy::too_many_arguments)]
fn map_query_to_reference_in_current_pool(
    reference_sketch: &ReferenceSketch,
    query_file: &QueryFile,
    kmer_size: usize,
    window_size: usize,
    min_identity: f64,
    mash_confidence: f64,
    frequency_threshold: usize,
    collect_metrics: bool,
) -> MappingOutput {
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
        .map(
            |(_scratch, mapping_results, mapping_metrics)| MappingOutput {
                results: mapping_results,
                metrics: mapping_metrics,
            },
        )
        .reduce(
            || MappingOutput {
                results: Vec::new(),
                metrics: MappingMetrics::default(),
            },
            |mut left, mut right| {
                left.results.append(&mut right.results);
                left.metrics.merge(right.metrics);
                left
            },
        )
}

/// Map a bounded batch of queries in one executor installation.
///
/// The outer parallel iterator lets small queries occupy workers together, while the nested
/// fragment iterator still spreads a single large query across the same run-scoped pool. Rayon
/// preserves the input query order in the collected result vector.
#[allow(clippy::too_many_arguments)]
pub(crate) fn map_query_batch_to_reference_parallel(
    executor: &MappingExecutor,
    reference_sketch: &ReferenceSketch,
    query_files: &[&QueryFile],
    kmer_size: usize,
    window_size: usize,
    min_identity: f64,
    mash_confidence: f64,
    frequency_threshold: usize,
    collect_metrics: bool,
) -> io::Result<Vec<MappingOutput>> {
    let Some(pool) = executor.pool.as_ref() else {
        return Ok(query_files
            .iter()
            .map(|query_file| {
                map_query_to_reference(
                    reference_sketch,
                    query_file,
                    kmer_size,
                    window_size,
                    min_identity,
                    mash_confidence,
                    frequency_threshold,
                    collect_metrics,
                )
            })
            .collect());
    };

    Ok(pool.install(|| {
        query_files
            .par_iter()
            .map(|query_file| {
                map_query_to_reference_in_current_pool(
                    reference_sketch,
                    query_file,
                    kmer_size,
                    window_size,
                    min_identity,
                    mash_confidence,
                    frequency_threshold,
                    collect_metrics,
                )
            })
            .collect()
    }))
}

#[cfg(debug_assertions)]
fn record_candidate_discovery_metrics(
    mapping_metrics: &mut MappingMetrics,
    seed_hit_count: usize,
    candidate_region_count: usize,
) {
    mapping_metrics.candidate_discovery_calls += 1;
    mapping_metrics.seed_hits_collected += seed_hit_count;
    let _ = candidate_region_count;
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
    let candidate_discovery_start: Instant = Instant::now();
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
    mapping_metrics.candidate_discovery_elapsed += candidate_discovery_start.elapsed();
    #[cfg(debug_assertions)]
    if collect_metrics {
        record_candidate_discovery_metrics(
            mapping_metrics,
            scratch.seed_hits.len(),
            scratch.candidate_regions.len(),
        );
    }
    let candidate_region_count: usize = scratch.candidate_regions.len();
    mapping_metrics.candidate_regions_found += candidate_region_count;
    mapping_metrics.candidate_regions_scored += candidate_region_count;

    scratch.fragment_mappings.clear();

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
pub(crate) fn compact_reciprocal_best_mappings(
    mut mapping_results: Vec<MappingResult>,
    fragment_length: u32,
) -> Vec<MappingResult> {
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

    summary_mappings
}

/// Collapse raw fragment mappings into final per-reference-file ANI summaries.
pub(crate) fn final_ani_computation(
    mapping_results: Vec<MappingResult>,
    query_file: &QueryFile,
    reference_file_count: usize,
    fragment_length: u32,
) -> AniComputation {
    let summary_mappings = compact_reciprocal_best_mappings(mapping_results, fragment_length);

    let mut summaries: Vec<AniSummary> = (0..reference_file_count)
        .map(|_| AniSummary::default())
        .collect::<Vec<_>>();
    let mut contig_summaries: Vec<Vec<ContigAniSummary>> = (0..reference_file_count)
        .map(|_| {
            (0..query_file.contig_names.len())
                .map(|_| ContigAniSummary::default())
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let mut fragment_identities: Vec<Vec<f64>> = (0..reference_file_count)
        .map(|_| Vec::new())
        .collect::<Vec<_>>();
    let mut contig_fragment_identities: Vec<Vec<Vec<f64>>> = (0..reference_file_count)
        .map(|_| {
            (0..query_file.contig_names.len())
                .map(|_| Vec::new())
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let mut reciprocal_best_keys: HashSet<MappingResultKey> = HashSet::new();

    for fragment in &query_file.fragments {
        for per_reference in &mut contig_summaries {
            if let Some(contig_summary) = per_reference.get_mut(fragment.contig_id) {
                contig_summary.eligible_fragments += 1;
            }
        }
    }

    for mapping in summary_mappings {
        reciprocal_best_keys.insert(MappingResultKey::from_mapping(&mapping));
        let summary = &mut summaries[mapping.reference_file_id];
        summary.shared_fragments += 1;
        summary.shared_bases += u64::from(mapping.query_fragment_length);
        summary.weighted_identity_sum +=
            mapping.identity * f64::from(mapping.query_fragment_length);
        fragment_identities[mapping.reference_file_id].push(mapping.identity);
        if let Some(query_fragment) = query_file.fragments.get(mapping.query_fragment_id) {
            let contig_summary =
                &mut contig_summaries[mapping.reference_file_id][query_fragment.contig_id].summary;
            contig_summary.shared_fragments += 1;
            contig_summary.shared_bases += u64::from(mapping.query_fragment_length);
            contig_summary.weighted_identity_sum +=
                mapping.identity * f64::from(mapping.query_fragment_length);
            contig_fragment_identities[mapping.reference_file_id][query_fragment.contig_id]
                .push(mapping.identity);
        }
    }

    for (summary, identities) in summaries.iter_mut().zip(&fragment_identities) {
        if !identities.is_empty() {
            summary.distribution_stats = compute_distribution_stats(identities);
        }
    }
    for (per_reference_summary, per_reference_identities) in
        contig_summaries.iter_mut().zip(&contig_fragment_identities)
    {
        for (contig_summary, identities) in per_reference_summary
            .iter_mut()
            .zip(per_reference_identities)
        {
            if !identities.is_empty() {
                contig_summary.summary.distribution_stats = compute_distribution_stats(identities);
            }
        }
    }

    AniComputation {
        summaries,
        contig_summaries,
        reciprocal_best_keys,
    }
}

fn compare_query_bucket(left: &MappingResult, right: &MappingResult) -> Ordering {
    left.reference_file_id
        .cmp(&right.reference_file_id)
        .then_with(|| left.query_fragment_id.cmp(&right.query_fragment_id))
        .then_with(|| left.identity.total_cmp(&right.identity))
        .then_with(|| left.reference_contig_id.cmp(&right.reference_contig_id))
        .then_with(|| left.reference_start.cmp(&right.reference_start))
        .then_with(|| left.query_fragment_length.cmp(&right.query_fragment_length))
        .then_with(|| left.query_minimizer_count.cmp(&right.query_minimizer_count))
        .then_with(|| {
            left.reference_minimizer_count
                .cmp(&right.reference_minimizer_count)
        })
        .then_with(|| left.shared_minimizers.cmp(&right.shared_minimizers))
        .then_with(|| left.union_minimizers.cmp(&right.union_minimizers))
        .then_with(|| left.jaccard.total_cmp(&right.jaccard))
}

fn compare_refbin_bucket(
    left: &MappingResult,
    right: &MappingResult,
    fragment_length: u32,
) -> Ordering {
    left.reference_contig_id
        .cmp(&right.reference_contig_id)
        .then_with(|| {
            reference_position_bin(left.reference_start, fragment_length).cmp(
                &reference_position_bin(right.reference_start, fragment_length),
            )
        })
        .then_with(|| left.identity.total_cmp(&right.identity))
        .then_with(|| left.reference_file_id.cmp(&right.reference_file_id))
        .then_with(|| left.query_fragment_id.cmp(&right.query_fragment_id))
        .then_with(|| left.reference_start.cmp(&right.reference_start))
        .then_with(|| left.query_fragment_length.cmp(&right.query_fragment_length))
        .then_with(|| left.query_minimizer_count.cmp(&right.query_minimizer_count))
        .then_with(|| {
            left.reference_minimizer_count
                .cmp(&right.reference_minimizer_count)
        })
        .then_with(|| left.shared_minimizers.cmp(&right.shared_minimizers))
        .then_with(|| left.union_minimizers.cmp(&right.union_minimizers))
        .then_with(|| left.jaccard.total_cmp(&right.jaccard))
}

fn reference_position_bin(position: u32, fragment_length: u32) -> u32 {
    position / fragment_length.saturating_sub(20).max(1)
}

pub(crate) fn compute_distribution_stats(fragment_identities: &[f64]) -> AniDistributionStats {
    let count: usize = fragment_identities.len();
    if count == 0 {
        return AniDistributionStats::default();
    }

    let mut sorted_identities: Vec<f64> = fragment_identities.to_vec();
    sorted_identities.sort_by(f64::total_cmp);

    let median: f64 = median_from_sorted(&sorted_identities);
    let mut absolute_deviations: Vec<f64> = fragment_identities
        .iter()
        .map(|identity| (identity - median).abs())
        .collect();
    absolute_deviations.sort_by(f64::total_cmp);
    let mad: f64 = median_from_sorted(&absolute_deviations);
    let f99: f64 = fraction_at_or_above(fragment_identities, 99.0);
    let f80: f64 = fraction_at_or_below(fragment_identities, 80.0);

    let mean: f64 = fragment_identities.iter().sum::<f64>() / count as f64;
    if count == 1 {
        return AniDistributionStats {
            median,
            stddev: f64::NAN,
            mad,
            ci_95_lower: mean,
            ci_95_upper: mean,
            f99,
            p99: f64::NAN,
            f80,
            p80: f64::NAN,
        };
    }

    let variance: f64 = fragment_identities
        .iter()
        .map(|identity| {
            let delta: f64 = identity - mean;
            delta * delta
        })
        .sum::<f64>()
        / (count - 1) as f64;
    let stddev: f64 = variance.sqrt();
    let standard_error: f64 = stddev / (count as f64).sqrt();
    let ci_delta: f64 = t_critical_95(count - 1) * standard_error;

    AniDistributionStats {
        median,
        stddev,
        mad,
        ci_95_lower: mean - ci_delta,
        ci_95_upper: mean + ci_delta,
        f99,
        p99: upper_tail_binomial_p_value(fragment_identities, mean, mad, 99.0),
        f80,
        p80: lower_tail_binomial_p_value(fragment_identities, mean, mad, 80.0),
    }
}

fn median_from_sorted(sorted_values: &[f64]) -> f64 {
    let count: usize = sorted_values.len();
    if count % 2 == 1 {
        sorted_values[count / 2]
    } else {
        f64::midpoint(sorted_values[(count / 2) - 1], sorted_values[count / 2])
    }
}

fn fraction_at_or_above(fragment_identities: &[f64], threshold: f64) -> f64 {
    fragment_identities
        .iter()
        .filter(|identity| **identity >= threshold)
        .count() as f64
        / fragment_identities.len() as f64
}

fn fraction_at_or_below(fragment_identities: &[f64], threshold: f64) -> f64 {
    fragment_identities
        .iter()
        .filter(|identity| **identity <= threshold)
        .count() as f64
        / fragment_identities.len() as f64
}

fn upper_tail_binomial_p_value(
    fragment_identities: &[f64],
    mean: f64,
    mad: f64,
    threshold: f64,
) -> f64 {
    if !mad.is_finite() || mad <= 0.0 || fragment_identities.len() < 2 {
        return f64::NAN;
    }

    let observed_at_or_above: usize = fragment_identities
        .iter()
        .filter(|identity| **identity >= threshold)
        .count();
    let count: usize = fragment_identities.len();
    let t_statistic: f64 = (threshold - mean) / mad;
    let tail_probability: f64 = 1.0 - t_cdf_approx(t_statistic, (count - 1) as f64);

    binomial_survival(observed_at_or_above, tail_probability, count)
}

fn lower_tail_binomial_p_value(
    fragment_identities: &[f64],
    mean: f64,
    mad: f64,
    threshold: f64,
) -> f64 {
    if !mad.is_finite() || mad <= 0.0 || fragment_identities.len() < 2 {
        return f64::NAN;
    }

    let observed_at_or_below: usize = fragment_identities
        .iter()
        .filter(|identity| **identity <= threshold)
        .count();
    let count: usize = fragment_identities.len();
    let t_statistic: f64 = (threshold - mean) / mad;
    let tail_probability: f64 = t_cdf_approx(t_statistic, (count - 1) as f64);

    binomial_survival(observed_at_or_below, tail_probability, count)
}

fn t_critical_95(degrees_of_freedom: usize) -> f64 {
    const T_CRITICAL_95: [f64; 30] = [
        12.706, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365, 2.306, 2.262, 2.228, 2.201, 2.179, 2.160,
        2.145, 2.131, 2.120, 2.110, 2.101, 2.093, 2.086, 2.080, 2.074, 2.069, 2.064, 2.060, 2.056,
        2.052, 2.048, 2.045, 2.042,
    ];

    T_CRITICAL_95
        .get(degrees_of_freedom.saturating_sub(1))
        .copied()
        .unwrap_or(1.96)
}

#[cfg(test)]
mod tests {
    use super::{
        checked_query_coordinate, compact_reciprocal_best_mappings, compute_distribution_stats,
        final_ani_computation,
    };
    use crate::ani::model::{AniSummary, MappingResult, QueryFile, QueryFragment};
    use std::io;

    fn assert_close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() < 1e-12,
            "expected {expected}, got {actual}"
        );
    }

    fn mapping(
        reference_file_id: usize,
        reference_contig_id: usize,
        query_fragment_id: usize,
        reference_start: u32,
        identity: f64,
    ) -> MappingResult {
        MappingResult {
            reference_file_id,
            reference_contig_id,
            query_fragment_id,
            query_fragment_length: 100,
            reference_start,
            identity,
            query_minimizer_count: 10,
            reference_minimizer_count: 11,
            shared_minimizers: 9,
            union_minimizers: 12,
            jaccard: 0.75,
        }
    }

    fn assert_summary_eq(actual: &AniSummary, expected: &AniSummary) {
        assert_eq!(actual.shared_fragments, expected.shared_fragments);
        assert_eq!(actual.shared_bases, expected.shared_bases);
        assert_eq!(
            actual.weighted_identity_sum.to_bits(),
            expected.weighted_identity_sum.to_bits()
        );
        let actual_stats = actual.distribution_stats;
        let expected_stats = expected.distribution_stats;
        assert_eq!(
            actual_stats.median.to_bits(),
            expected_stats.median.to_bits()
        );
        assert_eq!(
            actual_stats.stddev.to_bits(),
            expected_stats.stddev.to_bits()
        );
        assert_eq!(actual_stats.mad.to_bits(), expected_stats.mad.to_bits());
        assert_eq!(
            actual_stats.ci_95_lower.to_bits(),
            expected_stats.ci_95_lower.to_bits()
        );
        assert_eq!(
            actual_stats.ci_95_upper.to_bits(),
            expected_stats.ci_95_upper.to_bits()
        );
        assert_eq!(actual_stats.f99.to_bits(), expected_stats.f99.to_bits());
        assert_eq!(actual_stats.f80.to_bits(), expected_stats.f80.to_bits());
    }

    #[test]
    fn per_partition_compaction_preserves_final_summaries_and_keys() {
        let query = QueryFile {
            fragments: vec![
                QueryFragment {
                    id: 0,
                    contig_id: 0,
                    start: 0,
                    end: 100,
                    length: 100,
                    minimizers: vec![1],
                    seed_minimizers: vec![1],
                },
                QueryFragment {
                    id: 1,
                    contig_id: 1,
                    start: 0,
                    end: 100,
                    length: 100,
                    minimizers: vec![2],
                    seed_minimizers: vec![2],
                },
            ],
            contig_names: vec!["query-a".to_owned(), "query-b".to_owned()],
            mapped_length: 200,
        };
        let partition_one = vec![
            mapping(0, 0, 0, 100, 90.0),
            mapping(0, 0, 0, 500, 95.0),
            mapping(0, 0, 1, 501, 96.0),
        ];
        let partition_two = vec![
            mapping(1, 1, 0, 0, 92.0),
            mapping(1, 1, 0, 100, 91.0),
            mapping(1, 1, 1, 1, 93.0),
        ];
        let raw = partition_one
            .iter()
            .chain(&partition_two)
            .cloned()
            .collect();
        let per_partition_compacted = compact_reciprocal_best_mappings(partition_one, 100)
            .into_iter()
            .chain(compact_reciprocal_best_mappings(partition_two, 100))
            .collect::<Vec<_>>();
        assert_eq!(per_partition_compacted.len(), 2);

        let expected = final_ani_computation(raw, &query, 2, 100);
        let actual = final_ani_computation(per_partition_compacted, &query, 2, 100);

        assert_eq!(actual.summaries.len(), expected.summaries.len());
        for (actual, expected) in actual.summaries.iter().zip(&expected.summaries) {
            assert_summary_eq(actual, expected);
        }
        for (actual_reference, expected_reference) in actual
            .contig_summaries
            .iter()
            .zip(&expected.contig_summaries)
        {
            for (actual, expected) in actual_reference.iter().zip(expected_reference) {
                assert_eq!(actual.eligible_fragments, expected.eligible_fragments);
                assert_summary_eq(&actual.summary, &expected.summary);
            }
        }
        assert!(actual.reciprocal_best_keys == expected.reciprocal_best_keys);
    }

    #[test]
    fn reciprocal_best_compaction_breaks_exact_location_ties_deterministically() {
        let lower = mapping(0, 0, 0, 100, 95.0);
        let mut higher = lower.clone();
        higher.jaccard = 0.80;

        let forward = compact_reciprocal_best_mappings(vec![lower.clone(), higher.clone()], 100);
        let reverse = compact_reciprocal_best_mappings(vec![higher, lower], 100);

        assert_eq!(forward.len(), 1);
        assert_eq!(reverse.len(), 1);
        assert_eq!(forward[0].jaccard.to_bits(), 0.80_f64.to_bits());
        assert_eq!(forward[0].jaccard.to_bits(), reverse[0].jaccard.to_bits());
    }

    #[test]
    fn distribution_stats_compute_median_for_odd_count() {
        let stats = compute_distribution_stats(&[1.0, 2.0, 3.0, 4.0, 5.0]);

        assert_close(stats.median, 3.0);
    }

    #[test]
    fn distribution_stats_compute_median_for_even_count() {
        let stats = compute_distribution_stats(&[1.0, 2.0, 3.0, 4.0]);

        assert_close(stats.median, 2.5);
    }

    #[test]
    fn distribution_stats_compute_sample_stddev() {
        let stats = compute_distribution_stats(&[2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0]);

        assert_close(stats.stddev, 2.138_089_935_299_395);
    }

    #[test]
    fn distribution_stats_compute_median_absolute_deviation() {
        let stats = compute_distribution_stats(&[2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0]);

        assert_close(stats.mad, 0.5);
    }

    #[test]
    fn distribution_stats_flag_many_high_ani_fragments_using_mad() {
        let mut identities: Vec<f64> = vec![97.0; 400];
        identities.extend(vec![98.0; 400]);
        identities.extend(vec![99.0; 200]);

        let stats = compute_distribution_stats(&identities);

        assert_close(stats.mad, 1.0);
        assert!(
            stats.p99 < 0.05,
            "expected significant P99, got {}",
            stats.p99
        );
    }

    #[test]
    fn distribution_stats_flag_many_low_ani_fragments_using_mad() {
        let mut identities: Vec<f64> = vec![80.0; 200];
        identities.extend(vec![81.0; 400]);
        identities.extend(vec![82.0; 400]);

        let stats = compute_distribution_stats(&identities);

        assert_close(stats.mad, 1.0);
        assert!(
            stats.p80 < 0.05,
            "expected significant P80, got {}",
            stats.p80
        );
    }

    #[cfg(target_pointer_width = "64")]
    #[test]
    fn query_coordinate_overflow_is_contextual_invalid_data() {
        let value = usize::try_from(u64::from(u32::MAX) + 1).unwrap();
        let error = checked_query_coordinate(value, "fragment end", "contig-A").unwrap_err();

        assert_eq!(error.kind(), io::ErrorKind::InvalidData);
        assert!(error.to_string().contains("contig-A"));
        assert!(error.to_string().contains("fragment end"));
        assert!(error.to_string().contains(&value.to_string()));
    }
}
