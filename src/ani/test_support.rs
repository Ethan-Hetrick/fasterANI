//! Shared fixtures for the crate's unit tests.

use crate::ani::{
    ShardManifest, ShardManifestEntry, DEFAULT_FRAGMENT_LENGTH, DEFAULT_KMER_SIZE,
    DEFAULT_MAX_SHARD_MINIMIZERS, DEFAULT_MINIMIZER_HASH_SEED, DEFAULT_MIN_FRAGMENT_LENGTH,
    DEFAULT_SPLIT_N_RUN, DEFAULT_WINDOW_SIZE, SKETCH_DATABASE_SCHEMA_VERSION, SKETCH_VERSION,
};

pub(crate) fn repeated_acgt(len: usize) -> Vec<u8> {
    (0..len).map(|i| b"ACGT"[i % 4]).collect()
}

pub(crate) fn sample_shard_manifest() -> ShardManifest {
    ShardManifest {
        sketch_format_version: SKETCH_VERSION,
        database_schema_version: SKETCH_DATABASE_SCHEMA_VERSION,
        k: DEFAULT_KMER_SIZE,
        w: DEFAULT_WINDOW_SIZE,
        minimizer_hash_seed: DEFAULT_MINIMIZER_HASH_SEED,
        fragment_length: DEFAULT_FRAGMENT_LENGTH,
        min_fragment_length: DEFAULT_MIN_FRAGMENT_LENGTH,
        split_n_run: DEFAULT_SPLIT_N_RUN,
        max_shard_minimizers: DEFAULT_MAX_SHARD_MINIMIZERS,
        total_references: 2,
        total_reference_contigs: 2,
        total_mapped_reference_length: 6000,
        total_reference_minimizers: 20,
        total_shard_unique_minimizers: 18,
        build_unix_seconds: 1,
        generation_id: "test-generation".to_string(),
        global_frequency_filename: "database.test-generation.frequencies.bin".to_string(),
        global_frequency_file_bytes: 1,
        total_unique_minimizers: 17,
        build_args: vec!["fasterANI".to_string()],
        reference_list_checksum: 42,
        shards: vec![ShardManifestEntry {
            shard_index: 1,
            filename: "database.1.fasketch".to_string(),
            first_reference: 0,
            reference_count: 2,
            reference_contigs: 2,
            mapped_reference_length: 6000,
            reference_minimizers: 20,
            unique_minimizers: 18,
            file_bytes: 1,
        }],
    }
}
