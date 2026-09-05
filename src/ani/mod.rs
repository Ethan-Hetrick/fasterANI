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

pub use pipeline::{run, run_started_at};
