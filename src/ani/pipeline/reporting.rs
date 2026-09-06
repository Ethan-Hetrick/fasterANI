//! Result aggregation and TSV reporting.

use std::{
    collections::HashSet,
    io::{self, Write},
    time::Instant,
};

use crate::ani::{
    mapping::final_ani_computation,
    model::{
        query::{
            AniComputation, AniSummary, ContigAniSummary, MappingResult, MappingResultKey,
            QueryFile, QueryFragment,
        },
        reference::{ReferenceContigName, ReferenceFile},
    },
};

pub(super) struct PairSummaryStats {
    pub(crate) emitted_pairs: usize,
    pub(crate) ani_min: f64,
    pub(crate) ani_max: f64,
    pub(crate) af_min: f64,
    pub(crate) af_max: f64,
}

impl Default for PairSummaryStats {
    fn default() -> Self {
        Self {
            emitted_pairs: 0,
            ani_min: f64::NAN,
            ani_max: f64::NAN,
            af_min: f64::NAN,
            af_max: f64::NAN,
        }
    }
}

impl PairSummaryStats {
    fn record_pair(&mut self, ani: f64, aligned_fraction: f64) {
        if self.emitted_pairs == 0 {
            self.ani_min = ani;
            self.ani_max = ani;
            self.af_min = aligned_fraction;
            self.af_max = aligned_fraction;
        } else {
            self.ani_min = self.ani_min.min(ani);
            self.ani_max = self.ani_max.max(ani);
            self.af_min = self.af_min.min(aligned_fraction);
            self.af_max = self.af_max.max(aligned_fraction);
        }
        self.emitted_pairs += 1;
    }

    pub(super) fn merge(&mut self, other: Self) {
        if other.emitted_pairs == 0 {
            return;
        }
        if self.emitted_pairs == 0 {
            *self = other;
            return;
        }

        self.emitted_pairs += other.emitted_pairs;
        self.ani_min = self.ani_min.min(other.ani_min);
        self.ani_max = self.ani_max.max(other.ani_max);
        self.af_min = self.af_min.min(other.af_min);
        self.af_max = self.af_max.max(other.af_max);
    }
}

pub(super) struct QueryOutputStats {
    pub(crate) summary_elapsed: std::time::Duration,
    pub(crate) pair_stats: PairSummaryStats,
}

pub(super) fn write_results_header(output: &mut dyn Write, per_contig: bool) -> io::Result<()> {
    if per_contig {
        writeln!(
            output,
            "query_file\treference_file\tquery_contig\teligible_fragments\tshared_fragments\tshared_bases\tANI\tmedian_ANI\tstddev\tMAD\tci_95_upper\tci_95_lower\tF99\tF80"
        )
    } else {
        writeln!(
            output,
            "query_file\treference_file\tANI\tAF\ttotal_fragments\tmedian_ANI\tstddev\tMAD\tci_95_upper\tci_95_lower\tF99\tF80"
        )
    }
}

fn aggregate_values(
    summary: &AniSummary,
    query_mapped_length: u64,
    fragment_length: u32,
) -> (f64, f64, f64) {
    let total_fragment_equivalents: f64 = query_mapped_length as f64 / f64::from(fragment_length);
    let aligned_fraction: f64 = if total_fragment_equivalents > 0.0 {
        let shared_fragment_equivalents: f64 =
            summary.shared_bases as f64 / f64::from(fragment_length);
        shared_fragment_equivalents / total_fragment_equivalents
    } else {
        f64::NAN
    };
    let ani: f64 = if summary.shared_bases > 0 {
        summary.weighted_identity_sum / summary.shared_bases as f64
    } else {
        f64::NAN
    };

    (ani, aligned_fraction, total_fragment_equivalents)
}

fn contig_ani(summary: &AniSummary) -> f64 {
    if summary.shared_bases > 0 {
        summary.weighted_identity_sum / summary.shared_bases as f64
    } else {
        f64::NAN
    }
}

fn write_aggregate_summary_comments(
    reference_files: &[ReferenceFile],
    summaries: &[AniSummary],
    query_path: &str,
    query_mapped_length: u64,
    fragment_length: u32,
    output: &mut dyn Write,
) -> io::Result<()> {
    writeln!(
        output,
        "# aggregate_summary_header\tquery_file\treference_file\tANI\tAF\ttotal_fragments\tmedian_ANI\tstddev\tMAD\tci_95_upper\tci_95_lower\tF99\tF80"
    )?;

    for (reference_file, summary) in reference_files.iter().zip(summaries) {
        let (ani, aligned_fraction, total_fragment_equivalents) =
            aggregate_values(summary, query_mapped_length, fragment_length);
        let stats = summary.distribution_stats;
        writeln!(
            output,
            "# aggregate_summary\t{query_path}\t{}\t{ani:.3}\t{aligned_fraction:.3}\t{total_fragment_equivalents:.2}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}",
            reference_file.path,
            stats.median,
            stats.stddev,
            stats.mad,
            stats.ci_95_upper,
            stats.ci_95_lower,
            stats.f99,
            stats.f80,
        )?;
    }

    Ok(())
}

pub(super) fn write_mapping_stats_header(output: &mut dyn Write) -> io::Result<()> {
    writeln!(
        output,
        "query_file\treference_file\tquery_contig\treference_contig\tquery_fragment_id\tquery_start\tquery_end\treference_start\treference_end\tidentity\tquery_minimizer_count\treference_minimizer_count\tshared_minimizers\tunion_minimizers\tjaccard\tfragment_length\tis_reciprocal_best"
    )
}

pub(super) fn compare_mapping_stats_rows(
    left: &MappingResult,
    right: &MappingResult,
) -> std::cmp::Ordering {
    left.reference_file_id
        .cmp(&right.reference_file_id)
        .then_with(|| left.query_fragment_id.cmp(&right.query_fragment_id))
        .then_with(|| left.reference_contig_id.cmp(&right.reference_contig_id))
        .then_with(|| left.reference_start.cmp(&right.reference_start))
        .then_with(|| left.query_fragment_length.cmp(&right.query_fragment_length))
        .then_with(|| left.identity.total_cmp(&right.identity))
        .then_with(|| left.query_minimizer_count.cmp(&right.query_minimizer_count))
        .then_with(|| {
            left.reference_minimizer_count
                .cmp(&right.reference_minimizer_count)
        })
        .then_with(|| left.shared_minimizers.cmp(&right.shared_minimizers))
        .then_with(|| left.union_minimizers.cmp(&right.union_minimizers))
        .then_with(|| left.jaccard.total_cmp(&right.jaccard))
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
    let mut ordered_mappings: Vec<&MappingResult> = mapping_results.iter().collect();
    ordered_mappings.sort_by(|left, right| compare_mapping_stats_rows(left, right));

    for mapping in ordered_mappings {
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
            .map_or("unknown", String::as_str);

        let reference_contig_opt: Option<&str> = reference_contig_names
            .and_then(|contigs| contigs.get(mapping.reference_contig_id))
            .map(|c| c.name.as_str());

        let reference_offset: u64 = reference_contig_names
            .and_then(|contigs| contigs.get(mapping.reference_contig_id))
            .map_or(0, |contig| u64::from(contig.segment_start));
        let reference_start: u64 =
            reference_offset.saturating_add(u64::from(mapping.reference_start));
        let reference_end: u64 =
            reference_start.saturating_add(u64::from(mapping.query_fragment_length));
        let is_reciprocal_best: bool =
            reciprocal_best_keys.contains(&MappingResultKey::from_mapping(mapping));

        write!(
            output,
            "{query_path}\t{}\t{query_contig}\t",
            reference_file.path
        )?;

        if let Some(name) = reference_contig_opt {
            write!(output, "{name}")?;
        } else {
            write!(output, "{}", mapping.reference_contig_id)?;
        }

        writeln!(
            output,
            "\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{}\t{}\t{}\t{}\t{:.8}\t{}\t{}",
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

fn write_per_contig_results(
    reference_files: &[ReferenceFile],
    contig_summaries: &[Vec<ContigAniSummary>],
    query_file: &QueryFile,
    query_path: &str,
    output: &mut dyn Write,
) -> io::Result<()> {
    for (reference_file, per_contig) in reference_files.iter().zip(contig_summaries) {
        for (contig_id, contig_name) in query_file.contig_names.iter().enumerate() {
            let contig_summary = per_contig
                .get(contig_id)
                .cloned()
                .unwrap_or_else(ContigAniSummary::default);
            let summary = contig_summary.summary;
            let stats = summary.distribution_stats;
            writeln!(
                output,
                "{query_path}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}",
                reference_file.path,
                contig_name,
                contig_summary.eligible_fragments,
                summary.shared_fragments,
                summary.shared_bases,
                contig_ani(&summary),
                stats.median,
                stats.stddev,
                stats.mad,
                stats.ci_95_upper,
                stats.ci_95_lower,
                stats.f99,
                stats.f80,
            )?;
        }
    }

    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub(super) fn write_query_outputs(
    reference_files: &[ReferenceFile],
    reference_contig_names: Option<&[ReferenceContigName]>,
    mapping_results: Vec<MappingResult>,
    query_file: &QueryFile,
    query_path: &str,
    output: &mut dyn Write,
    mapping_stats_output: Option<&mut (dyn Write + '_)>,
    fragment_length: u32,
    per_contig: bool,
    emit_header: bool,
) -> io::Result<QueryOutputStats> {
    let summary_start: Instant = Instant::now();
    let raw_mapping_results: Option<Vec<MappingResult>> = mapping_stats_output
        .is_some()
        .then(|| mapping_results.clone());
    let ani_computation: AniComputation = final_ani_computation(
        mapping_results,
        query_file,
        reference_files.len(),
        fragment_length,
    );

    if let (Some(stats_output), Some(raw_mapping_results)) =
        (mapping_stats_output, raw_mapping_results.as_deref())
    {
        write_mapping_stats(
            reference_files,
            reference_contig_names,
            raw_mapping_results,
            &ani_computation.reciprocal_best_keys,
            query_file,
            query_path,
            stats_output,
        )?;
    }
    let summary_elapsed: std::time::Duration = summary_start.elapsed();

    let query_mapped_length: u64 = query_file.mapped_length();
    let mut pair_stats: PairSummaryStats = PairSummaryStats::default();

    if per_contig {
        write_aggregate_summary_comments(
            reference_files,
            &ani_computation.summaries,
            query_path,
            query_mapped_length,
            fragment_length,
            output,
        )?;
        if emit_header {
            write_results_header(output, per_contig)?;
        }
        write_per_contig_results(
            reference_files,
            &ani_computation.contig_summaries,
            query_file,
            query_path,
            output,
        )?;
    }

    if emit_header && !per_contig {
        write_results_header(output, per_contig)?;
    }

    for (reference_file, summary) in reference_files.iter().zip(ani_computation.summaries.iter()) {
        if summary.shared_fragments == 0 {
            continue;
        }

        let (ani, aligned_fraction, total_fragment_equivalents) =
            aggregate_values(summary, query_mapped_length, fragment_length);
        let stats = summary.distribution_stats;
        if !per_contig {
            writeln!(
                output,
                "{query_path}\t{}\t{ani:.3}\t{aligned_fraction:.3}\t{total_fragment_equivalents:.2}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}",
                reference_file.path,
                stats.median,
                stats.stddev,
                stats.mad,
                stats.ci_95_upper,
                stats.ci_95_lower,
                stats.f99,
                stats.f80,
            )?;
        }
        pair_stats.record_pair(ani, aligned_fraction);
    }

    Ok(QueryOutputStats {
        summary_elapsed,
        pair_stats,
    })
}
