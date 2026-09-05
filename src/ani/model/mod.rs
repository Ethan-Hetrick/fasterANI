//! Core data types describing references and queries.

mod query;
mod reference;

pub(crate) use query::{
    AniComputation, AniDistributionStats, AniSummary, ContigAniSummary, MappingResult,
    MappingResultKey, MappingScratch, QueryFile, QueryFragment, ReferenceCandidateRegion,
};
pub(crate) use reference::{
    CachedReferenceMetadata, ContigRecord, ReferenceContig, ReferenceContigName, ReferenceContigs,
    ReferenceFile, ReferenceIndex, ReferenceMemoryEstimate, ReferenceMinimizer, ReferenceSketch,
    SeedHit, ShardBuildResult, ShardManifest, ShardManifestEntry, ShardPlan, ShardedBuildOptions,
    SketchBuildStats, SketchParams, TransientReferenceIndex,
};
