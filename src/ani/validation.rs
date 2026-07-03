//! Validation and default-value helpers for CLI and runtime parameters.

use crate::ani::{AniError, DEFAULT_FRAGMENT_LENGTH, DEFAULT_MAX_SHARD_MINIMIZERS};

pub(crate) fn default_fragment_length() -> u32 {
    DEFAULT_FRAGMENT_LENGTH
}

pub(crate) fn default_max_shard_minimizers() -> usize {
    DEFAULT_MAX_SHARD_MINIMIZERS
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
