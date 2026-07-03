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
    database_build_parallelism, effective_index_build_mode, emit_progress,
    estimate_partitioned_shard_memory_bytes, legacy_sketch_path, manifest_path, memory_mib,
    plan_shards_by_minimizers, reference_list_checksum, shard_entry_path, shard_filename,
    shard_manifest_compatibility_error, shard_path, unix_timestamp_seconds,
    validate_max_shard_minimizers, FastaInput, IndexBuildMode, ReferenceSketch, RuntimeOptions,
    ShardBuildResult, ShardManifest, ShardManifestEntry, ShardPlan, ShardedBuildOptions,
    SketchBuildStats, SketchParams, SKETCH_DATABASE_SCHEMA_VERSION, SKETCH_KEY_MODE,
    SKETCH_VERSION,
};

/// Reference database opened by the CLI, either legacy single-sketch or manifest-backed shards.
pub(crate) enum SketchDatabase {
    Single(ReferenceSketch),
    Sharded {
        prefix: PathBuf,
        manifest: ShardManifest,
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
            tmp_dir,
            max_shard_minimizers,
            ..
        } = shard_opts;
        let Some(prefix) = sketch_prefix else {
            return Ok(Self::Single(ReferenceSketch::collect(
                references,
                params,
                runtime_options,
            )?));
        };

        validate_max_shard_minimizers(max_shard_minimizers)?;

        let manifest_path: PathBuf = manifest_path(prefix);
        if manifest_path.exists() {
            let manifest: ShardManifest = Self::load_manifest(prefix, params)?;
            return Ok(Self::Sharded {
                prefix: prefix.to_path_buf(),
                manifest,
            });
        }

        if let Some(legacy_path) = legacy_sketch_path(prefix) {
            return Ok(Self::Single(ReferenceSketch::load(
                &legacy_path,
                params,
                load_contig_names,
                tmp_dir,
                runtime_options,
            )?));
        }

        let manifest: ShardManifest =
            Self::build_sharded(references, params, prefix, shard_opts, runtime_options)?;

        Ok(Self::Sharded {
            prefix: prefix.to_path_buf(),
            manifest,
        })
    }

    pub(crate) fn build_sharded(
        references: &[FastaInput],
        params: SketchParams,
        prefix: &Path,
        shard_opts: ShardedBuildOptions<'_>,
        runtime_options: RuntimeOptions,
    ) -> io::Result<ShardManifest> {
        let SketchParams {
            kmer_size,
            window_size,
            fragment_length,
            min_fragment_length,
            split_n_run,
        } = params;
        let ShardedBuildOptions {
            tmp_dir,
            bgzip,
            max_shard_minimizers,
            index_build_mode,
            threads,
        } = shard_opts;
        validate_max_shard_minimizers(max_shard_minimizers)?;
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
                    "event=start\tmode=sharded\tprefix={}\treferences={}\tmax_shard_minimizers={max_shard_minimizers}\tthreads={}",
                    prefix.display(),
                    references.len(),
                    threads
                ),
                build_start,
            );
        }

        let shard_plans: Vec<ShardPlan> = plan_shards_by_minimizers(
            references,
            kmer_size,
            window_size,
            split_n_run,
            max_shard_minimizers,
            threads,
            runtime_options,
        )?;
        let build_parallelism: usize = database_build_parallelism(threads, &shard_plans);

        if runtime_options.progress_enabled {
            emit_progress(
                "database_build",
                &format!(
                    "event=shards_planned\tshards={}\tbuild_parallelism={build_parallelism}\tthreads={threads}\tmax_shard_minimizers={max_shard_minimizers}",
                    shard_plans.len()
                ),
                build_start,
            );
        }

        let completed_shards: Arc<AtomicUsize> = Arc::new(AtomicUsize::new(0));
        let shard_worker_threads: usize = threads.max(1).div_ceil(build_parallelism.max(1));
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(build_parallelism)
            .build()
            .map_err(|err| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("failed to initialize database build thread pool: {err}"),
                )
            })?;

        let mut shard_results: Vec<ShardBuildResult> = pool.install(|| {
            shard_plans
                .par_iter()
                .enumerate()
                .map(|(shard_offset, shard_plan)| {
                    let shard_index: usize = shard_offset + 1;
                    let first_reference: usize = shard_plan.first_reference;
                    let shard_end: usize = first_reference
                        .checked_add(shard_plan.reference_count)
                        .ok_or_else(|| {
                            io::Error::new(
                                io::ErrorKind::InvalidData,
                                "shard reference range overflow",
                            )
                        })?;
                    let reference_chunk: &[FastaInput] = &references[first_reference..shard_end];
                    let shard_path: PathBuf = shard_path(prefix, shard_index, bgzip);

                    if runtime_options.progress_enabled {
                        let effective_index_build_mode: IndexBuildMode =
                            effective_index_build_mode(
                                index_build_mode,
                                shard_plan.estimated_minimizers,
                            );
                        emit_progress(
                            "database_build",
                            &format!(
                                "event=shard_start\tshard={shard_index}\tfirst_reference={first_reference}\treference_count={}\testimated_minimizers={}\testimated_memory_mib={:.3}\tindex_build_mode={}\trequested_index_build_mode={}\tpath={}",
                                reference_chunk.len(),
                                shard_plan.estimated_minimizers,
                                memory_mib(estimate_partitioned_shard_memory_bytes(
                                    shard_plan.estimated_minimizers
                                )),
                                effective_index_build_mode.name(),
                                index_build_mode.name(),
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
                        bgzip,
                        shard_plan.estimated_minimizers,
                        index_build_mode,
                        runtime_options.with_worker_threads(shard_worker_threads),
                    )?;

                    let shards_done: usize =
                        completed_shards.fetch_add(1, AtomicOrdering::Relaxed) + 1;
                    if runtime_options.progress_enabled {
                        emit_progress(
                            "database_build",
                            &format!(
                                "event=shard_complete\tshard={shard_index}\tshards_done={shards_done}\tshards_total={}\treferences_done={}\treferences_total={}\tcontigs={}\treference_minimizers={}",
                                shard_plans.len(),
                                first_reference + stats.reference_count,
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
                            filename: shard_filename(prefix, shard_index, bgzip),
                            first_reference,
                            reference_count: stats.reference_count,
                            reference_contigs: stats.reference_contig_count,
                            mapped_reference_length: stats.mapped_reference_length,
                            reference_minimizers: stats.reference_minimizer_count,
                            unique_minimizers: stats.unique_minimizer_count,
                        },
                    })
                })
                .collect::<io::Result<Vec<_>>>()
        })?;

        shard_results.sort_by_key(|result| result.entry.shard_index);
        let shards: Vec<ShardManifestEntry> = shard_results
            .into_iter()
            .map(|result| result.entry)
            .collect();
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

        let manifest: ShardManifest = ShardManifest {
            sketch_format_version: SKETCH_VERSION,
            database_schema_version: SKETCH_DATABASE_SCHEMA_VERSION,
            k: kmer_size,
            w: window_size,
            key_mode: SKETCH_KEY_MODE.to_string(),
            fragment_length,
            min_fragment_length,
            split_n_run,
            dust_enabled: false,
            max_shard_minimizers,
            total_references: references.len(),
            total_reference_contigs,
            total_mapped_reference_length,
            total_reference_minimizers,
            total_shard_unique_minimizers,
            build_unix_seconds: unix_timestamp_seconds()?,
            build_args: env::args().collect(),
            reference_list_checksum: reference_list_checksum(references),
            shards,
        };

        Self::write_manifest(prefix, &manifest)?;
        if runtime_options.progress_enabled {
            emit_progress(
                "database_build",
                &format!(
                    "event=complete\tmanifest={}\tshards={}\treferences={}\tcontigs={}\treference_minimizers={}",
                    manifest_path(prefix).display(),
                    manifest.shards.len(),
                    manifest.total_references,
                    manifest.total_reference_contigs,
                    manifest.total_reference_minimizers
                ),
                build_start,
            );
        }

        Ok(manifest)
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
        fs::write(path, manifest_bytes)
    }

    pub(crate) fn load_manifest(prefix: &Path, params: SketchParams) -> io::Result<ShardManifest> {
        let SketchParams {
            kmer_size,
            window_size,
            fragment_length,
            min_fragment_length,
            split_n_run,
        } = params;
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

        if let Some(error) = shard_manifest_compatibility_error(
            &manifest,
            kmer_size,
            window_size,
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
            if !path.exists() {
                return Err(io::Error::new(
                    io::ErrorKind::NotFound,
                    format!(
                        "sharded sketch manifest references missing shard: {}",
                        path.display()
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
            Self::Sharded { manifest, .. } => manifest.total_shard_unique_minimizers,
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
        sample_shard_manifest, shard_manifest_compatibility_error, ShardManifest,
        DEFAULT_FRAGMENT_LENGTH, DEFAULT_KMER_SIZE, DEFAULT_MIN_FRAGMENT_LENGTH,
        DEFAULT_SPLIT_N_RUN, DEFAULT_WINDOW_SIZE,
    };
    use std::io;

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
    fn shard_manifest_rejects_legacy_dust_database() {
        let mut manifest: ShardManifest = sample_shard_manifest();
        manifest.dust_enabled = true;

        let error: Option<String> = shard_manifest_compatibility_error(
            &manifest,
            DEFAULT_KMER_SIZE,
            DEFAULT_WINDOW_SIZE,
            DEFAULT_FRAGMENT_LENGTH,
            DEFAULT_MIN_FRAGMENT_LENGTH,
            DEFAULT_SPLIT_N_RUN,
        );

        assert!(error
            .expect("dust-enabled manifest should be rejected")
            .contains("removed --dust filter"));
    }
}
