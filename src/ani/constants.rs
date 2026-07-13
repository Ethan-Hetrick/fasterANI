//! Tunable algorithm defaults and on-disk sketch-format constants shared across the crate.

use rustc_hash::FxHashMap;

use crate::ani::SeedHit;

pub(crate) const DEFAULT_KMER_SIZE: usize = 16;
pub(crate) const DEFAULT_WINDOW_SIZE: usize = 24;
pub(crate) const MINIMIZER_HASH_SEED: u32 = 42;
pub(crate) const DEFAULT_FRAGMENT_LENGTH: u32 = 3000;
pub(crate) const DEFAULT_FRAGMENT_STRIDE: u32 = DEFAULT_FRAGMENT_LENGTH;
pub(crate) const DEFAULT_MIN_FRAGMENT_LENGTH: u32 = DEFAULT_FRAGMENT_LENGTH;
/// Default for `--split-N`: minimum run-length of ambiguous `N` bases that forces
/// a contig to be split. `0` means "never split on N runs".
pub(crate) const DEFAULT_SPLIT_N_RUN: usize = 0;
pub(crate) const DEFAULT_MIN_PERCENT_IDENTITY: f64 = 80.0;
pub(crate) const DEFAULT_MASH_CONFIDENCE: f64 = 0.9;
pub(crate) const DEFAULT_FREQ_THRESHOLD_PERCENT: f64 = 0.0;
/// Default for `--mphf-gamma`: boomphf MPHF size/build-time tradeoff.
pub(crate) const DEFAULT_MPHF_GAMMA: f64 = 10.0;
/// Default for `--max-shard-minimizers`: per-shard reference minimizer budget.
/// 500M targets roughly 10 GiB `.fasketch` shards for GTDB-like bacterial genomes
/// using the observed storage model of about 21 bytes per reference minimizer.
pub(crate) const DEFAULT_MAX_SHARD_MINIMIZERS: usize = 500_000_000;
pub(crate) const PARTITIONED_INDEX_MINIMIZER_THRESHOLD: usize = 100_000_000;
pub(crate) const DEFAULT_PARTITION_TARGET_BYTES: usize = 512 * 1024 * 1024;
pub(crate) const MIN_PARTITION_COUNT: usize = 16;
pub(crate) const MAX_PARTITION_COUNT: usize = 4096;
pub(crate) const PARTITION_BUFFER_RECORDS: usize = 65_536;
pub(crate) const ESTIMATED_PARTITIONED_SHARD_BYTES_PER_MINIMIZER: usize = 16;
pub(crate) const SKETCH_MAGIC: &[u8; 8] = b"FANIIDX1";
pub(crate) const SKETCH_VERSION: u32 = 12;
pub(crate) const SKETCH_DATABASE_SCHEMA_VERSION: u32 = 2;
/// On-disk tag recording how minimizer keys are encoded: canonical (strand-agnostic)
/// k-mers, 2-bit-per-base packed, stored as `u32`. Persisted in the sketch manifest so
/// incompatible databases are rejected on load.
pub(crate) const SKETCH_KEY_MODE: &str = "canonical-2bit-u32";
pub(crate) const REFERENCE_PROGRESS_INTERVAL: usize = 1000;
pub(crate) const SKETCH_KEY_PACK_PROGRESS_INTERVAL: usize = 5_000_000;
#[cfg(test)]
pub(crate) const SKETCH_CONTIG_PACK_PROGRESS_INTERVAL: usize = 100_000;

pub(crate) type MinimizerKey = u32;
pub(crate) type ReferenceHitMap = FxHashMap<MinimizerKey, Vec<SeedHit>>;
