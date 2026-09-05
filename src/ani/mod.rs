//! Reference sketching, query mapping, and ANI summarization for `fasterANI`.
//!
//! The module keeps the FastANI-style data flow explicit:
//! reference minimizers are indexed once, query fragments are mapped through
//! seed-hit candidate regions and sliding-window identity scoring, and final output is summarized
//! per query/reference file pair.

mod cli;
mod constants;
mod error;
mod io_util;
mod mapping;
mod mash;
mod metrics;
mod minimizer;
mod mmap;
mod model;
mod params_file;
mod pipeline;
mod runtime;
mod sketch;
#[cfg(test)]
mod test_support;
mod validation;

pub(crate) use cli::{parse_cli_args, CliArgs};
#[cfg(test)]
pub(crate) use constants::SKETCH_CONTIG_PACK_PROGRESS_INTERVAL;
pub(crate) use constants::{
    MinimizerKey, ReferenceHitMap, DEFAULT_FRAGMENT_LENGTH, DEFAULT_FRAGMENT_STRIDE,
    DEFAULT_FREQ_THRESHOLD_PERCENT, DEFAULT_KMER_SIZE, DEFAULT_MASH_CONFIDENCE,
    DEFAULT_MAX_SHARD_MINIMIZERS, DEFAULT_MINIMIZER_HASH_SEED, DEFAULT_MIN_FRAGMENT_LENGTH,
    DEFAULT_MIN_PERCENT_IDENTITY, DEFAULT_MPHF_GAMMA, DEFAULT_PARTITION_TARGET_BYTES,
    DEFAULT_SPLIT_N_RUN, DEFAULT_WINDOW_SIZE, ESTIMATED_PARTITIONED_SHARD_BYTES_PER_MINIMIZER,
    MAX_PARTITION_COUNT, MIN_PARTITION_COUNT, PARTITIONED_INDEX_MINIMIZER_THRESHOLD,
    PARTITION_BUFFER_RECORDS, REFERENCE_PROGRESS_INTERVAL, SKETCH_DATABASE_SCHEMA_VERSION,
    SKETCH_KEY_PACK_PROGRESS_INTERVAL, SKETCH_MAGIC, SKETCH_VERSION,
};
pub(crate) use error::AniError;
pub(crate) use io_util::{
    align_up, append_path_suffix, checked_section_end, compress_file_to_bgzf,
    decompress_to_scratch, is_gzip_path, is_stdin_path, open_fasta_reader, sketch_reference_name,
    slice_as_bytes, slice_as_bytes_mut, write_padding, FastaInput, ScratchFile,
};
pub(crate) use mapping::{
    compact_reciprocal_best_mappings, final_ani_computation, lower_bound_minimizer_position,
    map_query_batch_to_reference_parallel, MappingExecutor,
};
pub(crate) use mash::{
    binomial_survival, estimate_relaxed_minimum_shared_minimizers, fastani_mash_distance,
    mash_distance_lower_bound, t_cdf_approx,
};
pub(crate) use metrics::{
    MappingMetrics, MappingOutput, SEED_HIT_HISTOGRAM_OVERFLOW_LABEL,
    SEED_HIT_HISTOGRAM_UPPER_BOUNDS,
};
#[cfg(test)]
pub(crate) use minimizer::{
    canonical_minimizer_observation, select_seed_minimizers, usable_minimizer_window_count,
    MinimizerObservation,
};
pub(crate) use minimizer::{
    canonical_minimizers_with_positions, canonical_minimizers_with_super_kmers,
    expected_minimizer_window_count, fastani_compatible_fragment_mode,
    is_no_usable_fragments_error, mapped_length_from_fragment_ranges, query_fragment_ranges,
    query_fragment_sketch, split_sequence_ranges, QueryFragmentSketch, SlidingSketchCounter,
};
pub(crate) use mmap::{MmapFile, MmapReferenceContigs, MmapReferenceIndex};
pub(crate) use model::{
    AniComputation, AniDistributionStats, AniSummary, CachedReferenceMetadata, ContigAniSummary,
    ContigRecord, MappingResult, MappingResultKey, MappingScratch, QueryFile, QueryFragment,
    ReferenceCandidateRegion, ReferenceContig, ReferenceContigName, ReferenceContigs,
    ReferenceFile, ReferenceIndex, ReferenceMemoryEstimate, ReferenceMinimizer, ReferenceSketch,
    SeedHit, ShardBuildResult, ShardManifest, ShardManifestEntry, ShardPlan, ShardedBuildOptions,
    SketchBuildStats, SketchParams, TransientReferenceIndex,
};
pub(crate) use params_file::{load_params_file, ParamsFileConfig};
pub(crate) use runtime::{
    emit_progress, emit_runtime_progress, memory_mib, peak_rss_kb, performance_metrics_enabled,
    reference_build_struct_bytes, RuntimeOptions,
};
pub(crate) use sketch::{
    build_generation_id, build_global_frequency_artifact, database_build_parallelism,
    effective_index_build_mode, estimate_partitioned_shard_memory_bytes,
    for_each_extracted_reference_segment, global_frequency_entry_path, global_frequency_filename,
    global_frequency_path, legacy_sketch_path, manifest_path, partition_build_plan,
    plan_shards_by_minimizers, reference_list_checksum, shard_entry_path, shard_filename,
    shard_manifest_compatibility_error, shard_path, sidecar_entry_path, unix_timestamp_seconds,
    write_bytes_atomically, write_name_sidecar, GlobalFrequencyArtifactStats, GlobalFrequencyIndex,
    GroupedKeyRecord, IndexBuildMode, NameSidecar, PartitionBuildPlan, PartitionGroupResult,
    PartitionHitRecord, PartitionWriters, ReferenceExtractionStats, SketchDatabase, SketchOutput,
};
#[cfg(test)]
pub(crate) use sketch::{
    estimate_reference_minimizer_windows, estimate_selected_minimizers_from_windows,
    partition_id_for_key, plan_shards_from_minimizer_counts,
};
#[cfg(test)]
pub(crate) use test_support::{repeated_acgt, sample_shard_manifest};
pub(crate) use validation::{
    default_fragment_length, default_max_shard_minimizers, validate_fragment_length,
    validate_kmer_size, validate_mash_confidence, validate_mash_threshold,
    validate_max_shard_minimizers, validate_mphf_gamma, validate_window_size,
};

pub use pipeline::{run, run_started_at};
