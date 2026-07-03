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

pub(crate) use cli::*;
pub(crate) use constants::*;
pub(crate) use error::*;
pub(crate) use io_util::*;
pub(crate) use mapping::*;
pub(crate) use mash::*;
pub(crate) use metrics::*;
pub(crate) use minimizer::*;
pub(crate) use mmap::*;
pub(crate) use model::*;
pub(crate) use params_file::*;
pub(crate) use runtime::*;
pub(crate) use sketch::*;
#[cfg(test)]
pub(crate) use test_support::*;
pub(crate) use validation::*;

pub use pipeline::run;
