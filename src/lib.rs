//! `fasterANI`: a fast, FastANI-style average nucleotide identity estimator.
//!
//! The crate is organized around the [`ani`] module, which holds the full
//! sketching, indexing, query-mapping, and ANI-summarization pipeline. The
//! binary entry point is the thin wrapper in `src/main.rs`; library consumers
//! drive the same pipeline through [`run`].

pub mod ani;

pub use ani::{run, run_started_at};
