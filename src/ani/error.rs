//! Typed errors for invalid CLI/runtime parameters.
//!
//! [`AniError`] gives parameter-validation failures named, self-describing
//! variants instead of stringly-typed [`io::Error`]s. It converts back into an
//! `io::Error` with `ErrorKind::InvalidInput` at API boundaries, so the
//! existing `io::Result` signatures and user-facing messages are unchanged.

use std::io;

/// A validation failure for a user-supplied parameter.
#[derive(Debug, thiserror::Error)]
pub(crate) enum AniError {
    #[error("--max-shard-minimizers must be at least 1")]
    MaxShardMinimizersTooSmall,
    #[error("--kmer-size must be between 1 and 16 for canonical-2bit-u32 keys")]
    KmerSizeOutOfRange,
    #[error("--window-size must be at least 1")]
    WindowSizeTooSmall,
    #[error("--fragment-length must be at least 1")]
    FragmentLengthTooSmall,
    #[error("--mash-threshold must be a finite value between 0 and 100")]
    MashThresholdOutOfRange,
    #[error("--mash-confidence must be a finite value in [0, 1)")]
    MashConfidenceOutOfRange,
}

impl From<AniError> for io::Error {
    fn from(error: AniError) -> Self {
        io::Error::new(io::ErrorKind::InvalidInput, error.to_string())
    }
}
