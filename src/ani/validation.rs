//! Validation and default-value helpers for CLI and runtime parameters.

use crate::ani::{
    AniError, DEFAULT_FRAGMENT_LENGTH, DEFAULT_MAX_SHARD_MINIMIZERS,
    ESTIMATED_PARTITIONED_SHARD_BYTES_PER_MINIMIZER,
};

pub(crate) fn default_fragment_length() -> u32 {
    DEFAULT_FRAGMENT_LENGTH
}

pub(crate) fn default_max_shard_minimizers() -> usize {
    DEFAULT_MAX_SHARD_MINIMIZERS
}

pub(crate) fn default_max_shard_minimizers_for_runtime(
    threads: usize,
    max_memory_bytes: Option<u64>,
) -> usize {
    let Some(max_memory_bytes) = max_memory_bytes else {
        return DEFAULT_MAX_SHARD_MINIMIZERS;
    };

    let Ok(max_memory_bytes) = usize::try_from(max_memory_bytes) else {
        return DEFAULT_MAX_SHARD_MINIMIZERS;
    };

    let active_jobs: usize = threads.max(1);
    let per_job_bytes: usize = max_memory_bytes / active_jobs;
    let memory_sized_minimizers: usize =
        per_job_bytes / ESTIMATED_PARTITIONED_SHARD_BYTES_PER_MINIMIZER;

    memory_sized_minimizers.max(DEFAULT_MAX_SHARD_MINIMIZERS)
}

pub(crate) fn validate_max_shard_minimizers(max_shard_minimizers: usize) -> Result<(), AniError> {
    if max_shard_minimizers == 0 {
        return Err(AniError::MaxShardMinimizersTooSmall);
    }

    Ok(())
}

pub(crate) fn validate_kmer_size(kmer_size: usize) -> Result<(), AniError> {
    if !(1..=16).contains(&kmer_size) {
        return Err(AniError::KmerSizeOutOfRange);
    }

    Ok(())
}

pub(crate) fn validate_window_size(window_size: usize) -> Result<(), AniError> {
    if window_size == 0 {
        return Err(AniError::WindowSizeTooSmall);
    }

    Ok(())
}

pub(crate) fn validate_fragment_length(fragment_length: u32) -> Result<(), AniError> {
    if fragment_length == 0 {
        return Err(AniError::FragmentLengthTooSmall);
    }

    Ok(())
}

pub(crate) fn validate_mash_threshold(mash_threshold: f64) -> Result<(), AniError> {
    if !mash_threshold.is_finite() || !(0.0..=100.0).contains(&mash_threshold) {
        return Err(AniError::MashThresholdOutOfRange);
    }

    Ok(())
}

pub(crate) fn validate_mash_confidence(mash_confidence: f64) -> Result<(), AniError> {
    if !mash_confidence.is_finite() || !(0.0..1.0).contains(&mash_confidence) {
        return Err(AniError::MashConfidenceOutOfRange);
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::ani::{default_max_shard_minimizers_for_runtime, DEFAULT_MAX_SHARD_MINIMIZERS};

    #[test]
    fn default_max_shard_minimizers_scales_with_memory_and_threads() {
        let one_hundred_gib: u64 = 100 * 1024 * 1024 * 1024;
        let max_shard_minimizers: usize =
            default_max_shard_minimizers_for_runtime(12, Some(one_hundred_gib));

        assert!(max_shard_minimizers > DEFAULT_MAX_SHARD_MINIMIZERS);
        assert_eq!(
            default_max_shard_minimizers_for_runtime(12, None),
            DEFAULT_MAX_SHARD_MINIMIZERS
        );
    }
}
