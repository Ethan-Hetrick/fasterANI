//! Reference sketch construction, on-disk format, sharding, and database loading.

mod build;
mod database;
mod extract;
mod frequency;
mod lookup;
mod partition;
mod persist;
mod serialize;
mod stream;

pub(crate) use database::SketchDatabase;
pub(crate) use extract::{for_each_extracted_reference_segment, ReferenceExtractionStats};
pub(crate) use frequency::{
    build_global_frequency_artifact, GlobalFrequencyArtifactStats, GlobalFrequencyIndex,
};
pub(crate) use partition::{
    database_build_parallelism, effective_index_build_mode,
    estimate_partitioned_shard_memory_bytes, partition_build_plan, plan_shards_by_minimizers,
    shard_manifest_compatibility_error, GroupedKeyRecord, IndexBuildMode, PartitionBuildPlan,
    PartitionGroupResult, PartitionHitRecord, PartitionWriters,
};
#[cfg(test)]
pub(crate) use partition::{
    estimate_reference_minimizer_windows, estimate_selected_minimizers_from_windows,
    partition_id_for_key, plan_shards_from_minimizer_counts,
};
pub(crate) use serialize::{
    build_generation_id, global_frequency_entry_path, global_frequency_filename,
    global_frequency_path, legacy_sketch_path, manifest_path, reference_list_checksum,
    shard_entry_path, shard_filename, shard_path, sidecar_entry_path, unix_timestamp_seconds,
    write_bytes_atomically, write_name_sidecar, NameSidecar, SketchOutput,
};
