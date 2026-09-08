//! The `SketchDatabase` (single or sharded) and its load/build orchestration.

use std::{
    env, fs, io,
    path::{Path, PathBuf},
    sync::{
        atomic::{AtomicUsize, Ordering as AtomicOrdering},
        Arc,
    },
    time::Instant,
};

use rayon::prelude::*;

use crate::ani::{
    constants::{SKETCH_DATABASE_SCHEMA_VERSION, SKETCH_VERSION},
    io_util::{sketch_reference_name, FastaInput},
    model::reference::{
        ReferenceSketch, ShardBuildResult, ShardManifest, ShardManifestEntry, ShardPlan,
        ShardedBuildOptions, SketchBuildStats, SketchParams,
    },
    runtime::{emit_progress, memory_mib, RuntimeOptions},
    sketch::{
        frequency::{
            build_global_frequency_artifact, GlobalFrequencyArtifactStats, GlobalFrequencyIndex,
        },
        partition::{
            database_build_parallelism, estimate_partitioned_shard_memory_bytes,
            plan_shards_by_size, shard_manifest_compatibility_error,
        },
        serialize::{
            build_generation_id, legacy_sketch_path, manifest_path, reference_list_checksum,
            shard_entry_path, shard_filename, shard_path, unix_timestamp_seconds,
            write_bytes_atomically,
        },
    },
    validation::validate_max_shard_size_bytes,
};

/// Reference database opened by the CLI, either legacy single-sketch or manifest-backed shards.
pub(crate) enum SketchDatabase {
    Single(ReferenceSketch),
    Sharded {
        prefix: PathBuf,
        manifest: ShardManifest,
        global_frequencies: Arc<GlobalFrequencyIndex>,
    },
}

impl SketchDatabase {
    pub(crate) fn collect_or_load(
        references: &[FastaInput],
        params: SketchParams,
        sketch_prefix: Option<&Path>,
        shard_opts: ShardedBuildOptions<'_>,
        load_contig_names: bool,
        runtime_options: RuntimeOptions,
    ) -> io::Result<Self> {
        let ShardedBuildOptions {
            max_shard_size_bytes,
            ..
        } = shard_opts;
        let Some(prefix) = sketch_prefix else {
            return Ok(Self::Single(ReferenceSketch::collect(
                references,
                params,
                runtime_options,
            )?));
        };

        validate_max_shard_size_bytes(max_shard_size_bytes)?;

        let manifest_path: PathBuf = manifest_path(prefix);
        if manifest_path.exists() {
            let manifest: ShardManifest = Self::load_manifest(prefix, params)?;
            if !references.is_empty() {
                let supplied_checksum: u64 = reference_list_checksum(references);
                // Older shards retain basenames, not their original source paths.
                // Updated generations therefore validate ordered stored identities;
                // fresh/legacy generations keep their original path-checksum rule.
                let matches = if let Some(ids) = &manifest.reference_identifiers {
                    *ids == references
                        .iter()
                        .map(|r| crate::ani::io_util::sketch_reference_name(&r.output_label))
                        .collect::<Vec<_>>()
                } else {
                    supplied_checksum == manifest.reference_list_checksum
                };
                if references.len() != manifest.total_references || !matches {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!(
                            "reference sketch database does not match the supplied reference list: supplied_references={} cached_references={} supplied_checksum={supplied_checksum} cached_checksum={}; use a different sketch prefix or omit reference inputs to load the existing database",
                            references.len(),
                            manifest.total_references,
                            manifest.reference_list_checksum
                        ),
                    ));
                }
            }
            let global_frequencies: Arc<GlobalFrequencyIndex> =
                Arc::new(GlobalFrequencyIndex::load(prefix, &manifest)?);
            return Ok(Self::Sharded {
                prefix: prefix.to_path_buf(),
                manifest,
                global_frequencies,
            });
        }

        if let Some(legacy_path) = legacy_sketch_path(prefix) {
            let sketch: ReferenceSketch =
                ReferenceSketch::load(&legacy_path, params, load_contig_names, runtime_options)?;
            if !references.is_empty()
                && (references.len() != sketch.files.len()
                    || references
                        .iter()
                        .zip(&sketch.files)
                        .any(|(reference, cached)| {
                            sketch_reference_name(&reference.output_label) != cached.path
                        }))
            {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "reference sketch does not match the supplied reference list; use a different sketch prefix or omit reference inputs to load the existing cache",
                ));
            }
            return Ok(Self::Single(sketch));
        }

        let manifest: ShardManifest =
            Self::build_sharded(references, params, prefix, shard_opts, runtime_options)?;
        let global_frequencies: Arc<GlobalFrequencyIndex> =
            Arc::new(GlobalFrequencyIndex::load(prefix, &manifest)?);

        Ok(Self::Sharded {
            prefix: prefix.to_path_buf(),
            manifest,
            global_frequencies,
        })
    }

    pub(crate) fn build_sharded(
        references: &[FastaInput],
        params: SketchParams,
        prefix: &Path,
        shard_opts: ShardedBuildOptions<'_>,
        runtime_options: RuntimeOptions,
    ) -> io::Result<ShardManifest> {
        let _lock = super::update::DatabaseWriteLock::acquire(prefix)?;
        if manifest_path(prefix).exists() || legacy_sketch_path(prefix).is_some() {
            return Err(io::Error::new(
                io::ErrorKind::AlreadyExists,
                "database already exists; use update",
            ));
        }
        super::update::addition_ids(references)?;
        let build_start = Instant::now();
        let SketchParams {
            kmer_size,
            window_size,
            minimizer_hash_seed,
            fragment_length,
            min_fragment_length,
            split_n_run,
        } = params;
        let ShardedBuildOptions {
            max_shard_size_bytes,
            tmp_dir,
            ..
        } = shard_opts;
        let generation_id = build_generation_id()?;
        let shards = Self::build_shard_files(
            references,
            params,
            prefix,
            shard_opts,
            runtime_options,
            &generation_id,
            1,
        )?;
        let total_reference_contigs: usize =
            shards.iter().map(|shard| shard.reference_contigs).sum();
        let total_mapped_reference_length: u64 = shards
            .iter()
            .map(|shard| shard.mapped_reference_length)
            .sum();
        let total_reference_minimizers: usize =
            shards.iter().map(|shard| shard.reference_minimizers).sum();
        let total_shard_unique_minimizers: usize =
            shards.iter().map(|shard| shard.unique_minimizers).sum();
        let global_frequency_stats: GlobalFrequencyArtifactStats = build_global_frequency_artifact(
            prefix,
            &generation_id,
            &shards,
            params,
            tmp_dir,
            runtime_options,
        )?;

        let manifest: ShardManifest = ShardManifest {
            sketch_format_version: SKETCH_VERSION,
            database_schema_version: SKETCH_DATABASE_SCHEMA_VERSION,
            k: kmer_size,
            w: window_size,
            minimizer_hash_seed,
            fragment_length,
            min_fragment_length,
            split_n_run,
            max_shard_size_bytes,
            total_references: references.len(),
            total_reference_contigs,
            total_mapped_reference_length,
            total_reference_minimizers,
            total_shard_unique_minimizers,
            build_unix_seconds: unix_timestamp_seconds()?,
            generation_id,
            global_frequency_filename: global_frequency_stats.filename,
            global_frequency_file_bytes: global_frequency_stats.file_bytes,
            total_unique_minimizers: global_frequency_stats.unique_minimizers,
            build_args: env::args().collect(),
            reference_list_checksum: reference_list_checksum(references),
            reference_identifiers: None,
            shards,
        };

        Self::write_manifest(prefix, &manifest)?;
        if runtime_options.progress_enabled {
            emit_progress(
                "database_build",
                &format!(
                    "event=complete\tgeneration_id={}\tmanifest={}\tshards={}\treferences={}\tcontigs={}\treference_minimizers={}\tunique_minimizers={}\tglobal_frequency_file_bytes={}",
                    manifest.generation_id,
                    manifest_path(prefix).display(),
                    manifest.shards.len(),
                    manifest.total_references,
                    manifest.total_reference_contigs,
                    manifest.total_reference_minimizers,
                    manifest.total_unique_minimizers,
                    manifest.global_frequency_file_bytes
                ),
                build_start,
            );
        }

        Ok(manifest)
    }

    pub(crate) fn build_shard_files(
        references: &[FastaInput],
        params: SketchParams,
        prefix: &Path,
        shard_opts: ShardedBuildOptions<'_>,
        runtime_options: RuntimeOptions,
        generation_id: &str,
        first_shard_index: usize,
    ) -> io::Result<Vec<ShardManifestEntry>> {
        let SketchParams {
            kmer_size,
            window_size,
            split_n_run,
            ..
        } = params;
        let ShardedBuildOptions {
            tmp_dir,
            max_shard_size_bytes,
            threads,
            ..
        } = shard_opts;
        validate_max_shard_size_bytes(max_shard_size_bytes)?;
        if let Some(parent) = prefix
            .parent()
            .filter(|parent| !parent.as_os_str().is_empty())
        {
            fs::create_dir_all(parent)?;
        }

        let build_start: Instant = Instant::now();
        if runtime_options.progress_enabled {
            emit_progress(
                "database_build",
                &format!(
                    "event=start\tmode=sharded\tgeneration_id={generation_id}\tprefix={}\treferences={}\tmax_shard_size_bytes={max_shard_size_bytes}\texecutor_threads={}",
                    prefix.display(),
                    references.len(),
                    threads
                ),
                build_start,
            );
        }

        let shard_plans: Vec<ShardPlan> = plan_shards_by_size(
            references,
            kmer_size,
            window_size,
            split_n_run,
            max_shard_size_bytes,
            threads,
            runtime_options,
        )?;
        let build_parallelism: usize = database_build_parallelism(threads, &shard_plans);

        if runtime_options.progress_enabled {
            emit_progress(
                "database_build",
                &format!(
                    "event=shards_planned\tgeneration_id={generation_id}\tshards={}\tbuild_parallelism={build_parallelism}\texecutor_threads={threads}\tmax_shard_size_bytes={max_shard_size_bytes}",
                    shard_plans.len()
                ),
                build_start,
            );
        }

        let completed_shards: Arc<AtomicUsize> = Arc::new(AtomicUsize::new(0));
        let completed_references: Arc<AtomicUsize> = Arc::new(AtomicUsize::new(0));
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads.max(1))
            .build()
            .map_err(|err| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("failed to initialize database build thread pool: {err}"),
                )
            })?;

        let build_shard = |shard_offset: usize,
                           shard_plan: &ShardPlan|
         -> io::Result<ShardBuildResult> {
            let shard_index = first_shard_index.checked_add(shard_offset).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidInput, "shard index overflow")
            })?;
            let first_reference: usize = shard_plan.first_reference;
            let shard_end: usize = first_reference
                .checked_add(shard_plan.reference_count)
                .ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidData, "shard reference range overflow")
                })?;
            let reference_chunk: &[FastaInput] = &references[first_reference..shard_end];
            let shard_path: PathBuf = shard_path(prefix, &generation_id, shard_index);
            let shard_runtime_options: RuntimeOptions = runtime_options
                .with_worker_threads(threads)
                .with_build_progress(&generation_id, shard_index)?;

            if runtime_options.progress_enabled {
                emit_progress(
                        "database_build",
                        &format!(
                            "event=shard_start\tgeneration_id={generation_id}\tshard={shard_index}\tfirst_reference={first_reference}\treference_count={}\testimated_minimizers={}\testimated_file_bytes={}\testimated_memory_mib={:.3}\tbuild_strategy=external_memory\texecutor_threads={threads}\tpath={}",
                            reference_chunk.len(),
                            shard_plan.estimated_minimizers,
                            shard_plan.estimated_file_bytes,
                            memory_mib(estimate_partitioned_shard_memory_bytes(
                                shard_plan.estimated_minimizers
                            )),
                            shard_path.display()
                        ),
                        build_start,
                    );
            }

            let stats: SketchBuildStats = ReferenceSketch::collect_and_save_streaming(
                reference_chunk,
                params,
                &shard_path,
                tmp_dir,
                shard_plan.estimated_minimizers,
                shard_runtime_options,
            )?;
            let file_bytes: u64 = fs::metadata(&shard_path)?.len();

            let shards_done: usize = completed_shards.fetch_add(1, AtomicOrdering::Relaxed) + 1;
            let references_done: usize = completed_references
                .fetch_add(stats.reference_count, AtomicOrdering::Relaxed)
                + stats.reference_count;
            if runtime_options.progress_enabled {
                emit_progress(
                        "database_build",
                        &format!(
                            "event=shard_complete\tgeneration_id={generation_id}\tshard={shard_index}\tshards_done={shards_done}\tshards_total={}\treferences_done={references_done}\treferences_total={}\tcontigs={}\treference_minimizers={}",
                            shard_plans.len(),
                            references.len(),
                            stats.reference_contig_count,
                            stats.reference_minimizer_count
                        ),
                        build_start,
                    );
            }
            Ok(ShardBuildResult {
                entry: ShardManifestEntry {
                    shard_index,
                    filename: shard_filename(prefix, &generation_id, shard_index),
                    first_reference,
                    reference_count: stats.reference_count,
                    reference_contigs: stats.reference_contig_count,
                    mapped_reference_length: stats.mapped_reference_length,
                    reference_minimizers: stats.reference_minimizer_count,
                    unique_minimizers: stats.unique_minimizer_count,
                    estimated_file_bytes: shard_plan.estimated_file_bytes,
                    file_bytes,
                },
            })
        };

        let next_shard: AtomicUsize = AtomicUsize::new(0);
        let mut shard_results: Vec<ShardBuildResult> = pool.install(|| {
            (0..build_parallelism)
                .into_par_iter()
                .map(|_| {
                    let mut worker_results: Vec<ShardBuildResult> = Vec::new();
                    loop {
                        let shard_offset: usize = next_shard.fetch_add(1, AtomicOrdering::Relaxed);
                        let Some(shard_plan) = shard_plans.get(shard_offset) else {
                            break;
                        };
                        worker_results.push(build_shard(shard_offset, shard_plan)?);
                    }
                    Ok(worker_results)
                })
                .collect::<io::Result<Vec<Vec<ShardBuildResult>>>>()
                .map(|worker_results| worker_results.into_iter().flatten().collect())
        })?;

        shard_results.sort_by_key(|result| result.entry.shard_index);
        let shards: Vec<ShardManifestEntry> = shard_results
            .into_iter()
            .map(|result| result.entry)
            .collect();
        Ok(shards)
    }

    pub(crate) fn write_manifest(prefix: &Path, manifest: &ShardManifest) -> io::Result<()> {
        let path: PathBuf = manifest_path(prefix);
        if let Some(parent) = path
            .parent()
            .filter(|parent| !parent.as_os_str().is_empty())
        {
            fs::create_dir_all(parent)?;
        }

        let manifest_bytes: Vec<u8> = serde_json::to_vec_pretty(manifest).map_err(|err| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("failed to encode sharded sketch manifest: {err}"),
            )
        })?;
        write_bytes_atomically(&path, &manifest_bytes)
    }

    pub(crate) fn read_manifest(prefix: &Path) -> io::Result<ShardManifest> {
        let path: PathBuf = manifest_path(prefix);
        let manifest_bytes: Vec<u8> = fs::read(&path)?;
        let manifest: ShardManifest = serde_json::from_slice(&manifest_bytes).map_err(|err| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "failed to read sharded sketch manifest {}: {err}",
                    path.display()
                ),
            )
        })?;

        Ok(manifest)
    }

    pub(crate) fn load_manifest(prefix: &Path, params: SketchParams) -> io::Result<ShardManifest> {
        let SketchParams {
            kmer_size,
            window_size,
            minimizer_hash_seed,
            fragment_length,
            min_fragment_length,
            split_n_run,
        } = params;
        let manifest = Self::read_manifest(prefix)?;

        if let Some(error) = shard_manifest_compatibility_error(
            &manifest,
            kmer_size,
            window_size,
            minimizer_hash_seed,
            fragment_length,
            min_fragment_length,
            split_n_run,
        ) {
            return Err(io::Error::new(io::ErrorKind::InvalidData, error));
        }

        let reference_count_sum: usize = manifest
            .shards
            .iter()
            .map(|shard| shard.reference_count)
            .sum();
        if reference_count_sum != manifest.total_references {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "sharded sketch manifest reference count mismatch: shards={reference_count_sum} total={}",
                    manifest.total_references
                ),
            ));
        }

        for shard in &manifest.shards {
            let path: PathBuf = shard_entry_path(prefix, shard);
            let metadata: fs::Metadata = fs::metadata(&path).map_err(|err| {
                io::Error::new(
                    err.kind(),
                    format!(
                        "sharded sketch manifest references unreadable shard {}: {err}",
                        path.display()
                    ),
                )
            })?;
            if metadata.len() != shard.file_bytes {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "sharded sketch file size mismatch for {}: manifest={} actual={}; rebuild the sketch database",
                        path.display(),
                        shard.file_bytes,
                        metadata.len()
                    ),
                ));
            }
        }
        Ok(manifest)
    }

    pub(crate) fn reference_count(&self) -> usize {
        match self {
            Self::Single(sketch) => sketch.files.len(),
            Self::Sharded { manifest, .. } => manifest.total_references,
        }
    }

    pub(crate) fn contig_count(&self) -> usize {
        match self {
            Self::Single(sketch) => sketch.contigs.len(),
            Self::Sharded { manifest, .. } => manifest.total_reference_contigs,
        }
    }

    pub(crate) fn unique_minimizer_count(&self) -> usize {
        match self {
            Self::Single(sketch) => sketch.index.len(),
            Self::Sharded { manifest, .. } => manifest.total_unique_minimizers,
        }
    }

    pub(crate) fn mode_name(
        &self,
        sketch_was_requested: bool,
        existing_database_loaded: bool,
    ) -> &'static str {
        match self {
            Self::Single(_) if !sketch_was_requested => "build-fasta",
            Self::Single(_) if existing_database_loaded => "load-sketch",
            Self::Single(_) => "build-fasta",
            Self::Sharded { .. } if existing_database_loaded => "load-sharded-sketch",
            Self::Sharded { .. } => "build-sharded-sketch",
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::ani::{
        constants::{
            DEFAULT_FRAGMENT_LENGTH, DEFAULT_KMER_SIZE, DEFAULT_MAX_SHARD_SIZE_BYTES,
            DEFAULT_MINIMIZER_HASH_SEED, DEFAULT_MIN_FRAGMENT_LENGTH, DEFAULT_SPLIT_N_RUN,
            DEFAULT_WINDOW_SIZE,
        },
        io_util::{append_path_suffix, FastaInput},
        model::reference::{ReferenceSketch, ShardManifest, ShardedBuildOptions, SketchParams},
        runtime::RuntimeOptions,
        sketch::{
            database::SketchDatabase,
            partition::shard_manifest_compatibility_error,
            serialize::{
                global_frequency_path, manifest_path, reference_list_checksum, shard_entry_path,
                shard_path, NameSidecar,
            },
        },
        test_support::sample_shard_manifest,
    };
    use std::{
        env, fs, io,
        path::PathBuf,
        time::{SystemTime, UNIX_EPOCH},
    };

    fn default_params() -> SketchParams {
        SketchParams {
            kmer_size: DEFAULT_KMER_SIZE,
            window_size: DEFAULT_WINDOW_SIZE,
            minimizer_hash_seed: DEFAULT_MINIMIZER_HASH_SEED,
            fragment_length: DEFAULT_FRAGMENT_LENGTH,
            min_fragment_length: DEFAULT_MIN_FRAGMENT_LENGTH,
            split_n_run: DEFAULT_SPLIT_N_RUN,
        }
    }

    fn unique_test_dir(label: &str) -> PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock after epoch")
            .as_nanos();
        env::temp_dir().join(format!("fasterani-{label}-{}-{nanos}", std::process::id()))
    }

    fn write_test_fasta(path: &std::path::Path, seed: u64) -> io::Result<()> {
        let mut state = seed;
        let mut sequence = String::with_capacity(6_000);
        for _ in 0..6_000 {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            sequence.push(match state & 3 {
                0 => 'A',
                1 => 'C',
                2 => 'G',
                _ => 'T',
            });
        }
        fs::write(path, format!(">sequence\n{sequence}\n"))
    }

    #[test]
    fn shard_manifest_round_trips_json() -> io::Result<()> {
        let manifest: ShardManifest = sample_shard_manifest();
        let encoded: Vec<u8> = serde_json::to_vec(&manifest).map_err(|err| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("failed to encode test manifest: {err}"),
            )
        })?;
        let decoded: ShardManifest = serde_json::from_slice(&encoded).map_err(|err| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("failed to decode test manifest: {err}"),
            )
        })?;

        assert_eq!(decoded.total_references, manifest.total_references);
        assert_eq!(decoded.shards[0].filename, "database.1.fasketch");

        Ok(())
    }

    #[test]
    fn shard_manifest_rejects_a_different_minimizer_seed() {
        let manifest: ShardManifest = sample_shard_manifest();
        let requested_seed = DEFAULT_MINIMIZER_HASH_SEED.wrapping_add(1);

        let error = shard_manifest_compatibility_error(
            &manifest,
            DEFAULT_KMER_SIZE,
            DEFAULT_WINDOW_SIZE,
            requested_seed,
            DEFAULT_FRAGMENT_LENGTH,
            DEFAULT_MIN_FRAGMENT_LENGTH,
            DEFAULT_SPLIT_N_RUN,
        )
        .expect("seed mismatch should be rejected");

        assert!(error.contains("minimizer_hash_seed"));
        assert!(error.contains("rebuild"));
    }

    #[test]
    fn existing_manifest_rejects_changed_supplied_references() -> io::Result<()> {
        let directory = unique_test_dir("stale-reference-list");
        fs::create_dir_all(&directory)?;
        let prefix = directory.join("database");
        let shard_path = directory.join("database.1.fasketch");
        fs::write(&shard_path, [0u8])?;

        let cached_references = vec![
            FastaInput::from_path("cached-a.fa".to_string()),
            FastaInput::from_path("cached-b.fa".to_string()),
        ];
        let mut manifest = sample_shard_manifest();
        manifest.reference_list_checksum = reference_list_checksum(&cached_references);
        SketchDatabase::write_manifest(&prefix, &manifest)?;

        let supplied_references = vec![
            FastaInput::from_path("cached-a.fa".to_string()),
            FastaInput::from_path("changed-b.fa".to_string()),
        ];
        let error = match SketchDatabase::collect_or_load(
            &supplied_references,
            default_params(),
            Some(&prefix),
            ShardedBuildOptions {
                tmp_dir: None,
                max_shard_size_bytes: DEFAULT_MAX_SHARD_SIZE_BYTES,
                threads: 1,
            },
            false,
            RuntimeOptions::default(),
        ) {
            Ok(_) => panic!("changed supplied references should be rejected"),
            Err(error) => error,
        };

        assert_eq!(error.kind(), io::ErrorKind::InvalidInput);
        assert!(error.to_string().contains("does not match"));
        fs::remove_dir_all(directory)?;
        Ok(())
    }

    #[test]
    fn manifest_load_rejects_a_changed_shard_size() -> io::Result<()> {
        let directory = unique_test_dir("shard-size-mismatch");
        fs::create_dir_all(&directory)?;
        let prefix = directory.join("database");
        fs::write(directory.join("database.1.fasketch"), [0u8, 1u8])?;
        SketchDatabase::write_manifest(&prefix, &sample_shard_manifest())?;

        let error = SketchDatabase::load_manifest(&prefix, default_params())
            .err()
            .expect("changed shard size should be rejected");
        assert_eq!(error.kind(), io::ErrorKind::InvalidData);
        assert!(error.to_string().contains("file size mismatch"));

        fs::remove_dir_all(directory)?;
        Ok(())
    }

    #[test]
    fn unpublished_rebuild_generation_does_not_replace_the_published_database() -> io::Result<()> {
        let directory = unique_test_dir("interrupted-sharded-rebuild");
        fs::create_dir_all(&directory)?;
        let reference_path = directory.join("reference.fna");
        write_test_fasta(&reference_path, 0x1234_5678_9abc_def0)?;
        let references = vec![FastaInput::from_path(
            reference_path.to_string_lossy().into_owned(),
        )];
        let prefix = directory.join("database");
        let params = default_params();
        let build_options = ShardedBuildOptions {
            tmp_dir: Some(&directory),
            max_shard_size_bytes: DEFAULT_MAX_SHARD_SIZE_BYTES,
            threads: 1,
        };

        let published_database = SketchDatabase::collect_or_load(
            &references,
            params,
            Some(&prefix),
            build_options,
            false,
            RuntimeOptions::default(),
        )?;
        let (published_generation, published_shard) = match &published_database {
            SketchDatabase::Sharded {
                prefix, manifest, ..
            } => (
                manifest.generation_id.clone(),
                shard_entry_path(prefix, &manifest.shards[0]),
            ),
            SketchDatabase::Single(_) => panic!("expected a sharded database"),
        };
        drop(published_database);
        let published_manifest_bytes = fs::read(manifest_path(&prefix))?;

        let incomplete_generation = "unpublished-interrupted-generation";
        fs::write(
            shard_path(&prefix, incomplete_generation, 1),
            b"incomplete shard",
        )?;
        fs::write(
            global_frequency_path(&prefix, incomplete_generation),
            b"incomplete global frequencies",
        )?;
        fs::write(
            append_path_suffix(&manifest_path(&prefix), ".tmp.interrupted"),
            b"{\"generation_id\":\"unpublished-interrupted-generation\"",
        )?;

        let loaded_database = SketchDatabase::collect_or_load(
            &[],
            params,
            Some(&prefix),
            build_options,
            false,
            RuntimeOptions::default(),
        )?;
        match &loaded_database {
            SketchDatabase::Sharded { manifest, .. } => {
                assert_eq!(manifest.generation_id, published_generation);
                assert_ne!(manifest.generation_id, incomplete_generation);
            }
            SketchDatabase::Single(_) => panic!("expected a sharded database"),
        }
        assert_eq!(fs::read(manifest_path(&prefix))?, published_manifest_bytes);

        let loaded_shard =
            ReferenceSketch::load(&published_shard, params, false, RuntimeOptions::default())?;
        assert_eq!(loaded_shard.files.len(), references.len());
        assert_eq!(loaded_shard.files[0].path, "reference.fna");
        let name_sidecar_path = fs::read_dir(&directory)?
            .filter_map(|entry| entry.ok().map(|entry| entry.path()))
            .find(|path| {
                path.file_name()
                    .and_then(|name| name.to_str())
                    .is_some_and(|name| name.starts_with("fasterani-names."))
            })
            .expect("sharded build name sidecar");
        let name_sidecar =
            NameSidecar::open(&name_sidecar_path, fs::metadata(&name_sidecar_path)?.len())?;
        assert_eq!(name_sidecar.genome_count(), 1);
        assert_eq!(name_sidecar.genome_name(0)?, "reference.fna");
        assert_eq!(name_sidecar.genome_length(0)?, 6_000);
        assert_eq!(name_sidecar.contig_name(0)?, "sequence");

        drop(loaded_shard);
        drop(loaded_database);
        drop(name_sidecar);
        fs::remove_dir_all(directory)?;
        Ok(())
    }
}
