//! Validation and default-value helpers for CLI and runtime parameters.

use crate::ani::{
    constants::{DEFAULT_FRAGMENT_LENGTH, DEFAULT_MAX_SHARD_MINIMIZERS},
    error::AniError,
};

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
    if !mash_confidence.is_finite() || !(0.0..=1.0).contains(&mash_confidence) {
        return Err(AniError::MashConfidenceOutOfRange);
    }

    Ok(())
}

pub(crate) fn validate_mphf_gamma(mphf_gamma: f64) -> Result<(), AniError> {
    if !mphf_gamma.is_finite() || mphf_gamma <= 1.01 {
        return Err(AniError::MphfGammaOutOfRange);
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{validate_mash_confidence, validate_mphf_gamma};

    #[test]
    fn mash_confidence_accepts_closed_unit_interval() {
        validate_mash_confidence(0.0).expect("0 should be valid");
        validate_mash_confidence(1.0).expect("1 should be valid");
    }

    #[test]
    fn mash_confidence_rejects_values_outside_unit_interval() {
        assert!(validate_mash_confidence(-0.1).is_err());
        assert!(validate_mash_confidence(1.1).is_err());
        assert!(validate_mash_confidence(f64::NAN).is_err());
    }

    #[test]
    fn mphf_gamma_accepts_values_above_boomphf_minimum() {
        validate_mphf_gamma(1.7).expect("default gamma should be valid");
    }

    #[test]
    fn mphf_gamma_rejects_invalid_values() {
        assert!(validate_mphf_gamma(1.01).is_err());
        assert!(validate_mphf_gamma(1.0).is_err());
        assert!(validate_mphf_gamma(f64::NAN).is_err());
        assert!(validate_mphf_gamma(f64::INFINITY).is_err());
    }
}
