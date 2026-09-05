//! Streaming and partitioned sketch construction that writes caches directly to disk.

use std::{
    fs, io,
    io::{BufReader, BufWriter, Read, Write},
    mem::{align_of, size_of},
    path::{Path, PathBuf},
    sync::atomic::{AtomicU32, AtomicUsize, Ordering as AtomicOrdering},
    time::Instant,
};

use boomphf::Mphf;
use rayon::prelude::*;

use crate::ani::{
    constants::{
        MinimizerKey, ReferenceHitMap, PARTITION_BUFFER_RECORDS, REFERENCE_PROGRESS_INTERVAL,
        SKETCH_KEY_PACK_PROGRESS_INTERVAL, SKETCH_MAGIC, SKETCH_VERSION,
    },
    io_util::{
        align_up, checked_section_end, slice_as_bytes, slice_as_bytes_mut, write_padding,
        FastaInput, ScratchFile,
    },
    model::reference::{
        CachedReferenceMetadata, ContigRecord, ReferenceContigName, ReferenceFile,
        ReferenceMinimizer, ReferenceSketch, SeedHit, SketchBuildStats, SketchParams,
    },
    runtime::{emit_runtime_progress, memory_mib, RuntimeOptions},
    sketch::{
        extract::{for_each_extracted_reference_segment, ReferenceExtractionStats},
        partition::{
            effective_index_build_mode, partition_build_plan, GroupedKeyRecord, IndexBuildMode,
            PartitionBuildPlan, PartitionGroupResult, PartitionHitRecord, PartitionWriters,
        },
        serialize::{write_name_sidecar, SketchOutput},
    },
};

const GROUPED_KEY_PACK_CHUNK: usize = 1_000_000;

impl ReferenceSketch {
    /// Build a reference sketch cache while streaming contig minimizers through scratch files.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn collect_and_save_streaming(
        references: &[FastaInput],
        params: SketchParams,
        cache_path: &Path,
        tmp_dir: Option<&Path>,
        estimated_minimizers: usize,
        index_build_mode: IndexBuildMode,
        runtime_options: RuntimeOptions,
    ) -> io::Result<SketchBuildStats> {
        if rayon::current_thread_index().is_none() {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(runtime_options.effective_worker_threads())
                .build()
                .map_err(|err| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("failed to initialize reference build thread pool: {err}"),
                    )
                })?;
            return pool.install(|| {
                Self::collect_and_save_streaming_in_current_pool(
                    references,
                    params,
                    cache_path,
                    tmp_dir,
                    estimated_minimizers,
                    index_build_mode,
                    runtime_options,
                )
            });
        }

        Self::collect_and_save_streaming_in_current_pool(
            references,
            params,
            cache_path,
            tmp_dir,
            estimated_minimizers,
            index_build_mode,
            runtime_options,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn collect_and_save_streaming_in_current_pool(
        references: &[FastaInput],
        params: SketchParams,
        cache_path: &Path,
        tmp_dir: Option<&Path>,
        estimated_minimizers: usize,
        index_build_mode: IndexBuildMode,
        runtime_options: RuntimeOptions,
    ) -> io::Result<SketchBuildStats> {
        let effective_mode: IndexBuildMode =
            effective_index_build_mode(index_build_mode, estimated_minimizers);
        match effective_mode {
            IndexBuildMode::Hash => Self::collect_and_save_streaming_hash(
                references,
                params,
                cache_path,
                tmp_dir,
                runtime_options,
            ),
            IndexBuildMode::Partitioned => Self::collect_and_save_streaming_partitioned(
                references,
                params,
                cache_path,
                tmp_dir,
                estimated_minimizers,
                runtime_options,
            ),
            IndexBuildMode::Auto => unreachable!("auto mode is resolved before sketch build"),
        }
    }

    /// Build a reference sketch cache using the in-memory minimizer hit map.
    pub(crate) fn collect_and_save_streaming_hash(
        references: &[FastaInput],
        params: SketchParams,
        cache_path: &Path,
        tmp_dir: Option<&Path>,
        runtime_options: RuntimeOptions,
    ) -> io::Result<SketchBuildStats> {
        let SketchParams { split_n_run, .. } = params;
        let build_start: Instant = Instant::now();
        let mut files: Vec<ReferenceFile> = Vec::new();
        let mut contig_records: Vec<ContigRecord> = Vec::new();
        let mut contig_names: Vec<ReferenceContigName> = Vec::new();
        let mut index: ReferenceHitMap = ReferenceHitMap::default();
        let (reference_minimizer_scratch, reference_minimizer_file): (ScratchFile, fs::File) =
            ScratchFile::create(tmp_dir, "reference-build-minimizers")?;
        let mut reference_minimizer_writer: BufWriter<fs::File> =
            BufWriter::new(reference_minimizer_file);
        let mut reference_minimizer_count: usize = 0usize;
        let mut total_seed_hits: usize = 0usize;
        files.reserve(references.len());
        contig_records.reserve(references.len());

        if runtime_options.progress_enabled {
            emit_runtime_progress(
                runtime_options,
                "reference_build",
                &format!(
                    "event=start\tmode=streaming\tfiles_total={}\tsplit_n_run={split_n_run}\tartifact={}\ttmp={}",
                    references.len(),
                    cache_path.display(),
                    reference_minimizer_scratch.path.display()
                ),
                build_start,
            );
        }
        for (file_id, reference) in references.iter().enumerate() {
            let file_id_u32: u32 = u32::try_from(file_id).map_err(|err| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "reference {} file id exceeds sketch cache limit: {err}",
                        reference.output_label
                    ),
                )
            })?;
            let extraction: ReferenceExtractionStats = for_each_extracted_reference_segment(
                reference,
                params,
                |segment| {
                    let reference_contig_id_u32: u32 =
                        u32::try_from(contig_records.len()).map_err(|err| {
                            io::Error::new(
                                io::ErrorKind::InvalidData,
                                format!(
                                    "reference {} record {:?}: contig id exceeds sketch cache limit: {err}",
                                    reference.output_label, segment.record_name
                                ),
                            )
                        })?;
                    index.reserve(segment.minimizers.len());

                    for minimizer in &segment.minimizers {
                        index.entry(minimizer.hash).or_default().push(SeedHit {
                            reference_contig_id: reference_contig_id_u32,
                            position: minimizer.position,
                        });
                    }

                    let minimizer_offset: u64 =
                        u64::try_from(reference_minimizer_count).map_err(|err| {
                            io::Error::new(
                                io::ErrorKind::InvalidData,
                                format!(
                                    "reference {} minimizer offset exceeds u64: {err}",
                                    reference.output_label
                                ),
                            )
                        })?;
                    let minimizer_count: u32 =
                        u32::try_from(segment.minimizers.len()).map_err(|err| {
                            io::Error::new(
                                io::ErrorKind::InvalidData,
                                format!(
                                    "reference {} record {:?}: contig has too many minimizers for sketch cache: {err}",
                                    reference.output_label, segment.record_name
                                ),
                            )
                        })?;

                    reference_minimizer_writer.write_all(slice_as_bytes(&segment.minimizers))?;
                    reference_minimizer_count = reference_minimizer_count
                        .checked_add(segment.minimizers.len())
                        .ok_or_else(|| {
                            io::Error::new(
                                io::ErrorKind::InvalidData,
                                format!(
                                    "reference {} minimizer count exceeds usize",
                                    reference.output_label
                                ),
                            )
                        })?;
                    total_seed_hits = total_seed_hits
                        .checked_add(segment.minimizers.len())
                        .ok_or_else(|| {
                            io::Error::new(
                                io::ErrorKind::InvalidData,
                                format!(
                                    "reference {} hit count exceeds usize",
                                    reference.output_label
                                ),
                            )
                        })?;
                    contig_records.push(ContigRecord {
                        minimizer_offset,
                        file_id: file_id_u32,
                        minimizer_count,
                    });
                    contig_names.push(ReferenceContigName {
                        file_id,
                        name: segment.record_name,
                        segment_start: segment.segment_start,
                        segment_end: segment.segment_end,
                    });
                    Ok(())
                },
            )?;

            files.push(ReferenceFile {
                path: reference.output_label.clone(),
                mapped_length: extraction.mapped_length,
                original_length: extraction.original_length,
            });

            let files_done: usize = file_id + 1;
            if runtime_options.progress_enabled
                && (files_done.is_multiple_of(REFERENCE_PROGRESS_INTERVAL)
                    || files_done == references.len())
            {
                emit_runtime_progress(
                    runtime_options,
                    "reference_build",
                    &format!(
                        "event=files\tmode=streaming\tfiles_done={files_done}\tfiles_total={}\tcontigs={}\treference_minimizers={reference_minimizer_count}\tunique_minimizers={}\tseed_hits={total_seed_hits}\tartifact={}",
                        references.len(),
                        contig_records.len(),
                        index.len(),
                        cache_path.display()
                    ),
                    build_start,
                );
            }
        }

        reference_minimizer_writer.flush()?;
        drop(reference_minimizer_writer);

        if runtime_options.progress_enabled {
            emit_runtime_progress(
                runtime_options,
                "reference_build",
                &format!(
                    "event=complete\tmode=streaming\tfiles_done={}\tcontigs={}\treference_minimizers={reference_minimizer_count}\tunique_minimizers={}\tseed_hits={total_seed_hits}\tartifact={}",
                    references.len(),
                    contig_records.len(),
                    index.len(),
                    cache_path.display()
                ),
                build_start,
            );
        }

        let stats: SketchBuildStats = SketchBuildStats {
            reference_count: files.len(),
            reference_contig_count: contig_records.len(),
            mapped_reference_length: files.iter().map(|file| file.mapped_length).sum(),
            reference_minimizer_count,
            unique_minimizer_count: index.len(),
        };

        Self::save_streamed_cache(
            cache_path,
            params,
            files,
            &index,
            contig_records,
            contig_names,
            reference_minimizer_count,
            &reference_minimizer_scratch,
            runtime_options,
        )?;

        Ok(stats)
    }

    /// Build a reference sketch cache using disk-partitioned minimizer hit records.
    pub(crate) fn collect_and_save_streaming_partitioned(
        references: &[FastaInput],
        params: SketchParams,
        cache_path: &Path,
        tmp_dir: Option<&Path>,
        estimated_minimizers: usize,
        runtime_options: RuntimeOptions,
    ) -> io::Result<SketchBuildStats> {
        let SketchParams { split_n_run, .. } = params;
        let build_start: Instant = Instant::now();
        let partition_plan: PartitionBuildPlan = partition_build_plan(estimated_minimizers);
        let mut files: Vec<ReferenceFile> = Vec::new();
        let mut contig_records: Vec<ContigRecord> = Vec::new();
        let mut contig_names: Vec<ReferenceContigName> = Vec::new();
        let mut partition_writers: PartitionWriters = PartitionWriters::new(
            partition_plan.partition_count,
            tmp_dir,
            PARTITION_BUFFER_RECORDS,
        )?;
        let (reference_minimizer_scratch, reference_minimizer_file): (ScratchFile, fs::File) =
            ScratchFile::create(tmp_dir, "reference-build-minimizers")?;
        let mut reference_minimizer_writer: BufWriter<fs::File> =
            BufWriter::new(reference_minimizer_file);
        let mut reference_minimizer_count: usize = 0usize;
        let mut total_seed_hits: usize = 0usize;
        files.reserve(references.len());
        contig_records.reserve(references.len());

        if runtime_options.progress_enabled {
            emit_runtime_progress(
                runtime_options,
                "reference_build",
                &format!(
                    "event=start\tmode=streaming\tindex_build_mode=partitioned\tfiles_total={}\tsplit_n_run={split_n_run}\tpartitions={}\testimated_minimizers={estimated_minimizers}\testimated_record_mib={:.3}\ttarget_partition_mib={:.3}\tartifact={}\ttmp={}",
                    references.len(),
                    partition_plan.partition_count,
                    memory_mib(partition_plan.estimated_record_bytes),
                    memory_mib(partition_plan.target_partition_bytes),
                    cache_path.display(),
                    reference_minimizer_scratch.path.display()
                ),
                build_start,
            );
        }
        for (file_id, reference) in references.iter().enumerate() {
            let file_id_u32: u32 = u32::try_from(file_id).map_err(|err| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "reference {} file id exceeds sketch cache limit: {err}",
                        reference.output_label
                    ),
                )
            })?;
            let extraction: ReferenceExtractionStats = for_each_extracted_reference_segment(
                reference,
                params,
                |segment| {
                    let reference_contig_id_u32: u32 =
                        u32::try_from(contig_records.len()).map_err(|err| {
                            io::Error::new(
                                io::ErrorKind::InvalidData,
                                format!(
                                    "reference {} record {:?}: contig id exceeds sketch cache limit: {err}",
                                    reference.output_label, segment.record_name
                                ),
                            )
                        })?;

                    // Build the record slice for this contig in one pass, then push as a batch.
                    // Avoids per-record function call overhead and enables SIMD partition ID computation.
                    for minimizer_chunk in segment.minimizers.chunks(PARTITION_BUFFER_RECORDS) {
                        let partition_hit_records: Vec<PartitionHitRecord> = minimizer_chunk
                            .iter()
                            .map(|minimizer| PartitionHitRecord {
                                key: minimizer.hash,
                                hit: SeedHit {
                                    reference_contig_id: reference_contig_id_u32,
                                    position: minimizer.position,
                                },
                            })
                            .collect();
                        partition_writers.push_batch(&partition_hit_records)?;
                    }

                    let minimizer_offset: u64 =
                        u64::try_from(reference_minimizer_count).map_err(|err| {
                            io::Error::new(
                                io::ErrorKind::InvalidData,
                                format!(
                                    "reference {} minimizer offset exceeds u64: {err}",
                                    reference.output_label
                                ),
                            )
                        })?;
                    let minimizer_count: u32 =
                        u32::try_from(segment.minimizers.len()).map_err(|err| {
                            io::Error::new(
                                io::ErrorKind::InvalidData,
                                format!(
                                    "reference {} record {:?}: contig has too many minimizers for sketch cache: {err}",
                                    reference.output_label, segment.record_name
                                ),
                            )
                        })?;

                    reference_minimizer_writer.write_all(slice_as_bytes(&segment.minimizers))?;
                    reference_minimizer_count = reference_minimizer_count
                        .checked_add(segment.minimizers.len())
                        .ok_or_else(|| {
                            io::Error::new(
                                io::ErrorKind::InvalidData,
                                format!(
                                    "reference {} minimizer count exceeds usize",
                                    reference.output_label
                                ),
                            )
                        })?;
                    total_seed_hits = total_seed_hits
                        .checked_add(segment.minimizers.len())
                        .ok_or_else(|| {
                            io::Error::new(
                                io::ErrorKind::InvalidData,
                                format!(
                                    "reference {} hit count exceeds usize",
                                    reference.output_label
                                ),
                            )
                        })?;
                    contig_records.push(ContigRecord {
                        minimizer_offset,
                        file_id: file_id_u32,
                        minimizer_count,
                    });
                    contig_names.push(ReferenceContigName {
                        file_id,
                        name: segment.record_name,
                        segment_start: segment.segment_start,
                        segment_end: segment.segment_end,
                    });
                    Ok(())
                },
            )?;

            files.push(ReferenceFile {
                path: reference.output_label.clone(),
                mapped_length: extraction.mapped_length,
                original_length: extraction.original_length,
            });

            let files_done: usize = file_id + 1;
            if runtime_options.progress_enabled
                && (files_done.is_multiple_of(REFERENCE_PROGRESS_INTERVAL)
                    || files_done == references.len())
            {
                emit_runtime_progress(
                    runtime_options,
                    "reference_build",
                    &format!(
                        "event=files\tmode=streaming\tindex_build_mode=partitioned\tfiles_done={files_done}\tfiles_total={}\tcontigs={}\treference_minimizers={reference_minimizer_count}\tseed_hits={total_seed_hits}\tpartitions={}\tartifact={}",
                        references.len(),
                        contig_records.len(),
                        partition_plan.partition_count,
                        cache_path.display()
                    ),
                    build_start,
                );
            }
        }

        reference_minimizer_writer.flush()?;
        drop(reference_minimizer_writer);
        partition_writers.flush_all()?;

        if runtime_options.progress_enabled {
            emit_runtime_progress(
                runtime_options,
                "reference_build",
                &format!(
                    "event=complete\tmode=streaming\tindex_build_mode=partitioned\tfiles_done={}\tcontigs={}\treference_minimizers={reference_minimizer_count}\tseed_hits={total_seed_hits}\tpartitions={}\tartifact={}",
                    references.len(),
                    contig_records.len(),
                    partition_plan.partition_count,
                    cache_path.display()
                ),
                build_start,
            );
        }

        let reference_count: usize = files.len();
        let reference_contig_count: usize = contig_records.len();
        let mapped_reference_length: u64 = files.iter().map(|file| file.mapped_length).sum();
        let unique_minimizer_count: usize = Self::save_partitioned_streamed_cache(
            cache_path,
            params,
            files,
            contig_records,
            contig_names,
            reference_minimizer_count,
            &reference_minimizer_scratch,
            &partition_writers,
            partition_plan,
            tmp_dir,
            runtime_options,
        )?;

        Ok(SketchBuildStats {
            reference_count,
            reference_contig_count,
            mapped_reference_length,
            reference_minimizer_count,
            unique_minimizer_count,
        })
    }

    pub(crate) fn read_partition_records(path: &Path) -> io::Result<Vec<PartitionHitRecord>> {
        let byte_len: u64 = fs::metadata(path)?.len();
        let record_size: u64 = size_of::<PartitionHitRecord>() as u64;
        if !byte_len.is_multiple_of(record_size) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "partition file {} has {byte_len} bytes, not a multiple of {record_size}",
                    path.display()
                ),
            ));
        }

        let record_count: usize = usize::try_from(byte_len / record_size).map_err(|err| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("partition file has too many records to load: {err}"),
            )
        })?;
        let mut records: Vec<PartitionHitRecord> =
            vec![PartitionHitRecord::default(); record_count];
        if record_count > 0 {
            let mut reader: BufReader<fs::File> = BufReader::new(fs::File::open(path)?);
            reader.read_exact(slice_as_bytes_mut(&mut records))?;
        }

        Ok(records)
    }

    pub(crate) fn group_partition_records(
        partition_index: usize,
        partition_path: &Path,
        tmp_dir: Option<&Path>,
    ) -> io::Result<PartitionGroupResult> {
        let mut records: Vec<PartitionHitRecord> = Self::read_partition_records(partition_path)?;
        records.sort_unstable_by_key(|record| {
            (
                record.key,
                record.hit.reference_contig_id,
                record.hit.position,
            )
        });

        let (grouped_key_scratch, grouped_key_file): (ScratchFile, fs::File) = ScratchFile::create(
            tmp_dir,
            &format!("partition-{partition_index}-grouped-keys"),
        )?;
        let (hit_payload_scratch, hit_payload_file): (ScratchFile, fs::File) = ScratchFile::create(
            tmp_dir,
            &format!("partition-{partition_index}-hit-payloads"),
        )?;
        let mut grouped_key_writer: BufWriter<fs::File> = BufWriter::new(grouped_key_file);
        let mut hit_payload_writer: BufWriter<fs::File> = BufWriter::new(hit_payload_file);
        let mut keys: Vec<MinimizerKey> = Vec::with_capacity(records.len());
        let mut local_hit_count: usize = 0usize;
        let mut hit_buffer: Vec<SeedHit> = Vec::with_capacity(records.len().min(1_048_576));

        let mut group_start: usize = 0usize;
        while group_start < records.len() {
            let key: MinimizerKey = records[group_start].key;
            let mut group_end: usize = group_start + 1;
            while group_end < records.len() && records[group_end].key == key {
                group_end += 1;
            }

            let hit_count: usize = group_end - group_start;
            let hit_count_u32: u32 = u32::try_from(hit_count).map_err(|err| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("minimizer hit count exceeds sketch cache limit: {err}"),
                )
            })?;
            let hit_offset_u64: u64 = u64::try_from(local_hit_count).map_err(|err| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("minimizer hit offset exceeds sketch cache limit: {err}"),
                )
            })?;
            let grouped_record: GroupedKeyRecord = GroupedKeyRecord {
                key,
                hit_offset: hit_offset_u64,
                hit_count: hit_count_u32,
            };
            grouped_key_writer.write_all(slice_as_bytes(std::slice::from_ref(&grouped_record)))?;
            keys.push(key);

            for record in &records[group_start..group_end] {
                hit_buffer.push(record.hit);
                if hit_buffer.len() >= 1_048_576 {
                    hit_payload_writer.write_all(slice_as_bytes(&hit_buffer))?;
                    hit_buffer.clear();
                }
            }

            local_hit_count = local_hit_count.checked_add(hit_count).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "partition hit count overflow")
            })?;
            group_start = group_end;
        }

        if !hit_buffer.is_empty() {
            hit_payload_writer.write_all(slice_as_bytes(&hit_buffer))?;
        }
        grouped_key_writer.flush()?;
        hit_payload_writer.flush()?;
        drop(grouped_key_writer);
        drop(hit_payload_writer);

        let expected_grouped_key_bytes: u64 = (keys.len() * size_of::<GroupedKeyRecord>()) as u64;
        let actual_grouped_key_bytes: u64 = fs::metadata(&grouped_key_scratch.path)?.len();
        if actual_grouped_key_bytes != expected_grouped_key_bytes {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "partition grouped key scratch file has {actual_grouped_key_bytes} bytes, expected {expected_grouped_key_bytes}"
                ),
            ));
        }

        let expected_hit_payload_bytes: u64 = (local_hit_count * size_of::<SeedHit>()) as u64;
        let actual_hit_payload_bytes: u64 = fs::metadata(&hit_payload_scratch.path)?.len();
        if actual_hit_payload_bytes != expected_hit_payload_bytes {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "partition hit payload scratch file has {actual_hit_payload_bytes} bytes, expected {expected_hit_payload_bytes}"
                ),
            ));
        }

        Ok(PartitionGroupResult {
            partition_index,
            grouped_key_scratch,
            hit_payload_scratch,
            key_count: keys.len(),
            keys,
            hit_count: local_hit_count,
        })
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn save_partitioned_streamed_cache(
        path: &Path,
        params: SketchParams,
        files: Vec<ReferenceFile>,
        contig_records: Vec<ContigRecord>,
        contig_names: Vec<ReferenceContigName>,
        reference_minimizer_count: usize,
        reference_minimizer_scratch: &ScratchFile,
        partition_writers: &PartitionWriters,
        partition_plan: PartitionBuildPlan,
        tmp_dir: Option<&Path>,
        runtime_options: RuntimeOptions,
    ) -> io::Result<usize> {
        let SketchParams {
            kmer_size,
            window_size,
            minimizer_hash_seed,
            fragment_length,
            min_fragment_length,
            split_n_run,
        } = params;
        if let Some(parent) = path
            .parent()
            .filter(|parent| !parent.as_os_str().is_empty())
        {
            fs::create_dir_all(parent)?;
        }

        let save_start: Instant = Instant::now();
        if runtime_options.progress_enabled {
            emit_runtime_progress(
                runtime_options,
                "sketch_save",
                &format!(
                    "event=start\tmode=streaming\tindex_build_mode=partitioned\tpartitions={}\testimated_minimizers={}\ttarget_partition_mib={:.3}\tcontigs={}\tfiles={}\tartifact={}\ttmp={}",
                    partition_plan.partition_count,
                    partition_plan.estimated_record_bytes / size_of::<PartitionHitRecord>(),
                    memory_mib(partition_plan.target_partition_bytes),
                    contig_records.len(),
                    files.len(),
                    path.display(),
                    reference_minimizer_scratch.path.display()
                ),
                save_start,
            );
        }
        let partition_paths: Vec<PathBuf> = (0..partition_writers.partition_count())
            .map(|partition_index| partition_writers.path(partition_index).to_path_buf())
            .collect();
        let partition_parallelism: usize = runtime_options
            .effective_worker_threads()
            .min(partition_paths.len().max(1));
        let completed_partitions: AtomicUsize = AtomicUsize::new(0);
        let mut partition_results: Vec<PartitionGroupResult> = partition_paths
            .par_iter()
            .enumerate()
            .map(|(partition_index, partition_path)| {
                let result: PartitionGroupResult =
                    Self::group_partition_records(partition_index, partition_path, tmp_dir)?;
                let partitions_done: usize =
                    completed_partitions.fetch_add(1, AtomicOrdering::Relaxed) + 1;
                if runtime_options.progress_enabled
                    && (partitions_done.is_multiple_of(16)
                        || partitions_done == partition_paths.len())
                {
                    emit_runtime_progress(
                        runtime_options,
                        "sketch_save",
                        &format!(
                            "event=sort_group_partitions\tmode=streaming\tindex_build_mode=partitioned\tpartitions_done={partitions_done}\tpartitions={}\tpartition_parallelism={partition_parallelism}\tartifact={}",
                            partition_paths.len(),
                            path.display()
                        ),
                        save_start,
                    );
                }
                Ok(result)
            })
            .collect::<io::Result<Vec<_>>>()?;
        partition_results.sort_by_key(|result| result.partition_index);

        let total_unique_minimizers: usize = partition_results
            .iter()
            .map(|result| result.key_count)
            .sum();
        let mut keys: Vec<MinimizerKey> = Vec::with_capacity(total_unique_minimizers);
        let mut total_hits: usize = 0usize;
        for result in &mut partition_results {
            keys.extend(std::mem::take(&mut result.keys));
            total_hits = total_hits.checked_add(result.hit_count).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "total hit count overflow")
            })?;
        }

        if runtime_options.progress_enabled {
            emit_runtime_progress(
                runtime_options,
                "sketch_save",
                &format!(
                    "event=keys_collected\tmode=streaming\tindex_build_mode=partitioned\tkey_count={}\thit_count={total_hits}\tpartition_parallelism={partition_parallelism}\tartifact={}",
                    keys.len(),
                    path.display(),
                ),
                save_start,
            );
        }
        let key_count: usize = keys.len();
        let mphf: Mphf<MinimizerKey> = Mphf::new_parallel(runtime_options.mphf_gamma, &keys, None);
        drop(keys);
        if runtime_options.progress_enabled {
            emit_runtime_progress(
                runtime_options,
                "sketch_save",
                &format!(
                    "event=mphf_built\tmode=streaming\tindex_build_mode=partitioned\tkey_count={key_count}\tartifact={}",
                    path.display(),
                ),
                save_start,
            );
        }
        let slot_keys: Vec<AtomicU32> = (0..key_count).map(|_| AtomicU32::new(0)).collect();
        let hit_offsets: Vec<AtomicU32> = (0..key_count).map(|_| AtomicU32::new(0)).collect();
        let hit_counts: Vec<AtomicU32> = (0..key_count).map(|_| AtomicU32::new(0)).collect();
        let mut grouped_chunk: Vec<GroupedKeyRecord> = Vec::with_capacity(GROUPED_KEY_PACK_CHUNK);
        let mut grouped_records_done: usize = 0usize;
        let mut partition_hit_offset: u64 = 0u64;
        for result in &partition_results {
            let mut grouped_key_reader: BufReader<fs::File> =
                BufReader::new(fs::File::open(&result.grouped_key_scratch.path)?);
            let mut partition_records_done: usize = 0usize;
            while partition_records_done < result.key_count {
                let records_to_read: usize =
                    (result.key_count - partition_records_done).min(GROUPED_KEY_PACK_CHUNK);
                grouped_chunk.resize(records_to_read, GroupedKeyRecord::default());
                grouped_key_reader.read_exact(slice_as_bytes_mut(&mut grouped_chunk))?;

                grouped_chunk.par_iter().try_for_each(|grouped_record| {
                    let slot: usize = mphf.hash(&grouped_record.key) as usize;
                    let global_offset: u64 = partition_hit_offset
                        .checked_add(grouped_record.hit_offset)
                        .ok_or_else(|| {
                            io::Error::new(io::ErrorKind::InvalidData, "global hit offset overflow")
                        })?;
                    let global_offset_u32: u32 = u32::try_from(global_offset).map_err(|_| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!(
                                "shard has more than u32::MAX hit-payload entries ({global_offset}); \
                                 reduce --max-shard-minimizers below 4 294 967 295"
                            ),
                        )
                    })?;
                    // MPHF slots are unique for the collected key set, so these independent stores
                    // are race-free. Relaxed ordering is sufficient because the Rayon join below
                    // happens before the atomics are converted into ordinary arrays.
                    slot_keys[slot].store(grouped_record.key, AtomicOrdering::Relaxed);
                    hit_offsets[slot].store(global_offset_u32, AtomicOrdering::Relaxed);
                    hit_counts[slot].store(grouped_record.hit_count, AtomicOrdering::Relaxed);
                    Ok::<(), io::Error>(())
                })?;

                partition_records_done += records_to_read;
                grouped_records_done += records_to_read;
                if runtime_options.progress_enabled
                    && (grouped_records_done.is_multiple_of(SKETCH_KEY_PACK_PROGRESS_INTERVAL)
                        || grouped_records_done == key_count)
                {
                    emit_runtime_progress(
                        runtime_options,
                        "sketch_save",
                        &format!(
                            "event=pack_index\tmode=streaming\tindex_build_mode=partitioned\tkeys_done={grouped_records_done}\tkey_count={key_count}\tartifact={}",
                            path.display(),
                        ),
                        save_start,
                    );
                }
            }
            partition_hit_offset = partition_hit_offset
                .checked_add(u64::try_from(result.hit_count).map_err(|err| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("partition hit count exceeds u64: {err}"),
                    )
                })?)
                .ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidData, "partition hit offset overflow")
                })?;
        }
        let slot_keys: Vec<MinimizerKey> =
            slot_keys.into_iter().map(AtomicU32::into_inner).collect();
        let hit_offsets: Vec<u32> = hit_offsets.into_iter().map(AtomicU32::into_inner).collect();
        let hit_counts: Vec<u32> = hit_counts.into_iter().map(AtomicU32::into_inner).collect();

        let expected_reference_minimizer_bytes: u64 =
            (reference_minimizer_count * size_of::<ReferenceMinimizer>()) as u64;
        let actual_reference_minimizer_bytes: u64 =
            fs::metadata(&reference_minimizer_scratch.path)?.len();
        if actual_reference_minimizer_bytes != expected_reference_minimizer_bytes {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "reference minimizer scratch file has {actual_reference_minimizer_bytes} bytes, expected {expected_reference_minimizer_bytes}"
                ),
            ));
        }
        if contig_names.len() != contig_records.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "contig sidecar would have {} records, expected {}",
                    contig_names.len(),
                    contig_records.len()
                ),
            ));
        }
        let (name_sidecar_filename, name_sidecar_file_bytes) =
            write_name_sidecar(path, &files, &contig_names)?;
        let metadata: CachedReferenceMetadata = CachedReferenceMetadata {
            version: SKETCH_VERSION,
            k: kmer_size,
            w: window_size,
            minimizer_hash_seed,
            fragment_length,
            min_fragment_length,
            split_n_run,
            reference_count: files.len(),
            mphf,
            key_count: slot_keys.len(),
            hit_count: total_hits,
            contig_count: contig_records.len(),
            reference_minimizer_count,
            name_sidecar_filename,
            name_sidecar_file_bytes,
        };
        let metadata_bytes: Vec<u8> = serde_json::to_vec(&metadata).map_err(|err| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("failed to encode reference sketch metadata: {err}"),
            )
        })?;

        let metadata_end: usize = SKETCH_MAGIC
            .len()
            .checked_add(size_of::<u64>())
            .and_then(|offset| offset.checked_add(metadata_bytes.len()))
            .ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "metadata length overflow")
            })?;
        let slot_keys_offset: usize = align_up(metadata_end, 8);
        let hit_offsets_offset: usize = align_up(
            checked_section_end(slot_keys_offset, slot_keys.len(), size_of::<MinimizerKey>())?,
            align_of::<u32>(),
        );
        let hit_counts_offset: usize =
            checked_section_end(hit_offsets_offset, hit_offsets.len(), size_of::<u32>())?;
        let hit_payloads_offset: usize = align_up(
            checked_section_end(hit_counts_offset, hit_counts.len(), size_of::<u32>())?,
            align_of::<SeedHit>(),
        );
        let contig_records_offset: usize = align_up(
            checked_section_end(hit_payloads_offset, total_hits, size_of::<SeedHit>())?,
            align_of::<ContigRecord>(),
        );
        let reference_minimizers_offset: usize = align_up(
            checked_section_end(
                contig_records_offset,
                contig_records.len(),
                size_of::<ContigRecord>(),
            )?,
            align_of::<ReferenceMinimizer>(),
        );

        let mut sketch_output: SketchOutput = SketchOutput::create(path)?;
        let mut writer: &mut BufWriter<fs::File> = sketch_output.writer_mut()?;
        writer.write_all(SKETCH_MAGIC)?;
        writer.write_all(&(metadata_bytes.len() as u64).to_le_bytes())?;
        writer.write_all(&metadata_bytes)?;
        write_padding(&mut writer, slot_keys_offset - metadata_end)?;
        writer.write_all(slice_as_bytes(&slot_keys))?;
        write_padding(
            &mut writer,
            hit_offsets_offset
                - checked_section_end(
                    slot_keys_offset,
                    slot_keys.len(),
                    size_of::<MinimizerKey>(),
                )?,
        )?;
        writer.write_all(slice_as_bytes(&hit_offsets))?;
        writer.write_all(slice_as_bytes(&hit_counts))?;
        write_padding(
            &mut writer,
            hit_payloads_offset
                - checked_section_end(hit_counts_offset, hit_counts.len(), size_of::<u32>())?,
        )?;
        for result in &partition_results {
            let mut hit_payload_reader: fs::File =
                fs::File::open(&result.hit_payload_scratch.path)?;
            io::copy(&mut hit_payload_reader, &mut writer)?;
        }
        write_padding(
            &mut writer,
            contig_records_offset
                - checked_section_end(hit_payloads_offset, total_hits, size_of::<SeedHit>())?,
        )?;
        writer.write_all(slice_as_bytes(&contig_records))?;
        write_padding(
            &mut writer,
            reference_minimizers_offset
                - checked_section_end(
                    contig_records_offset,
                    contig_records.len(),
                    size_of::<ContigRecord>(),
                )?,
        )?;
        let mut reference_minimizer_reader: fs::File =
            fs::File::open(&reference_minimizer_scratch.path)?;
        io::copy(&mut reference_minimizer_reader, &mut writer)?;
        let output_file_bytes: u64 = sketch_output.finish()?;

        if runtime_options.progress_enabled {
            emit_runtime_progress(
                runtime_options,
                "sketch_save",
                &format!(
                    "event=complete\tmode=streaming\tindex_build_mode=partitioned\tpath={}\tfile_bytes={}",
                    path.display(),
                    output_file_bytes
                ),
                save_start,
            );
        }
        Ok(key_count)
    }
}

#[cfg(test)]
mod tests {
    use crate::ani::{
        constants::{
            DEFAULT_FRAGMENT_LENGTH, DEFAULT_KMER_SIZE, DEFAULT_MINIMIZER_HASH_SEED,
            DEFAULT_MIN_FRAGMENT_LENGTH, DEFAULT_SPLIT_N_RUN, DEFAULT_WINDOW_SIZE,
        },
        io_util::FastaInput,
        model::reference::{ReferenceMinimizer, ReferenceSketch, SketchBuildStats, SketchParams},
        runtime::RuntimeOptions,
        sketch::{partition::IndexBuildMode, serialize::NameSidecar},
    };
    use std::{
        env, fs, io,
        path::PathBuf,
        time::{SystemTime, UNIX_EPOCH},
    };

    fn assert_partitioned_matches_hash(
        sequence: &[u8],
        expected_key_parity: usize,
        label: &str,
    ) -> io::Result<()> {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock after epoch")
            .as_nanos();
        let directory: PathBuf = env::temp_dir().join(format!(
            "fasterani-partitioned-{label}-{}-{nanos}",
            std::process::id()
        ));
        fs::create_dir_all(&directory)?;
        let reference_path: PathBuf = directory.join("reference.fa");
        let hash_sketch_path: PathBuf = directory.join("hash.fasketch");
        let partitioned_sketch_path: PathBuf = directory.join("partitioned.fasketch");
        fs::write(
            &reference_path,
            format!(">ref\n{}\n", String::from_utf8_lossy(sequence)),
        )?;
        let references: Vec<FastaInput> = vec![FastaInput::from_path(
            reference_path.to_string_lossy().into_owned(),
        )];
        let build_runtime: RuntimeOptions = RuntimeOptions::default().with_worker_threads(4);

        let hash_stats: SketchBuildStats = ReferenceSketch::collect_and_save_streaming(
            &references,
            SketchParams {
                kmer_size: DEFAULT_KMER_SIZE,
                window_size: DEFAULT_WINDOW_SIZE,
                minimizer_hash_seed: DEFAULT_MINIMIZER_HASH_SEED,
                fragment_length: DEFAULT_FRAGMENT_LENGTH,
                min_fragment_length: DEFAULT_MIN_FRAGMENT_LENGTH,
                split_n_run: DEFAULT_SPLIT_N_RUN,
            },
            &hash_sketch_path,
            None,
            1,
            IndexBuildMode::Hash,
            build_runtime,
        )?;
        let partitioned_stats: SketchBuildStats = ReferenceSketch::collect_and_save_streaming(
            &references,
            SketchParams {
                kmer_size: DEFAULT_KMER_SIZE,
                window_size: DEFAULT_WINDOW_SIZE,
                minimizer_hash_seed: DEFAULT_MINIMIZER_HASH_SEED,
                fragment_length: DEFAULT_FRAGMENT_LENGTH,
                min_fragment_length: DEFAULT_MIN_FRAGMENT_LENGTH,
                split_n_run: DEFAULT_SPLIT_N_RUN,
            },
            &partitioned_sketch_path,
            None,
            hash_stats.reference_minimizer_count,
            IndexBuildMode::Partitioned,
            build_runtime,
        )?;
        assert_eq!(
            hash_stats.unique_minimizer_count % 2,
            expected_key_parity,
            "{label} fixture did not exercise the intended key-count parity"
        );
        let hash_sketch: ReferenceSketch = ReferenceSketch::load(
            &hash_sketch_path,
            SketchParams {
                kmer_size: DEFAULT_KMER_SIZE,
                window_size: DEFAULT_WINDOW_SIZE,
                minimizer_hash_seed: DEFAULT_MINIMIZER_HASH_SEED,
                fragment_length: DEFAULT_FRAGMENT_LENGTH,
                min_fragment_length: DEFAULT_MIN_FRAGMENT_LENGTH,
                split_n_run: DEFAULT_SPLIT_N_RUN,
            },
            false,
            RuntimeOptions::default(),
        )?;
        let partitioned_sketch: ReferenceSketch = ReferenceSketch::load(
            &partitioned_sketch_path,
            SketchParams {
                kmer_size: DEFAULT_KMER_SIZE,
                window_size: DEFAULT_WINDOW_SIZE,
                minimizer_hash_seed: DEFAULT_MINIMIZER_HASH_SEED,
                fragment_length: DEFAULT_FRAGMENT_LENGTH,
                min_fragment_length: DEFAULT_MIN_FRAGMENT_LENGTH,
                split_n_run: DEFAULT_SPLIT_N_RUN,
            },
            false,
            RuntimeOptions::default(),
        )?;

        assert_eq!(
            hash_stats.reference_minimizer_count,
            partitioned_stats.reference_minimizer_count
        );
        assert_eq!(
            hash_stats.unique_minimizer_count,
            partitioned_stats.unique_minimizer_count
        );
        assert_eq!(hash_sketch.files.len(), partitioned_sketch.files.len());
        assert_eq!(hash_sketch.contigs.len(), partitioned_sketch.contigs.len());
        assert_eq!(hash_sketch.files[0].original_length, sequence.len() as u64);
        assert_eq!(
            partitioned_sketch.files[0].original_length,
            sequence.len() as u64
        );

        let name_sidecars: Vec<PathBuf> = fs::read_dir(&directory)?
            .filter_map(|entry| entry.ok().map(|entry| entry.path()))
            .filter(|path| {
                path.file_name()
                    .and_then(|name| name.to_str())
                    .is_some_and(|name| name.starts_with("fasterani-names."))
            })
            .collect();
        assert_eq!(name_sidecars.len(), 1);
        let name_sidecar_bytes = fs::metadata(&name_sidecars[0])?.len();
        let name_sidecar = NameSidecar::open(&name_sidecars[0], name_sidecar_bytes)?;
        assert_eq!(name_sidecar.genome_count(), 1);
        assert_eq!(name_sidecar.genome_name(0)?, "reference.fa");
        assert_eq!(name_sidecar.genome_length(0)?, sequence.len() as u64);
        assert_eq!(name_sidecar.contig_count(), hash_sketch.contigs.len());
        assert_eq!(name_sidecar.contig_name(0)?, "ref");

        for contig_id in 0..hash_sketch.contigs.len() {
            let hash_minimizers: &[ReferenceMinimizer] = hash_sketch
                .contigs
                .minimizers(contig_id)
                .expect("hash contig");
            let partitioned_minimizers: &[ReferenceMinimizer] = partitioned_sketch
                .contigs
                .minimizers(contig_id)
                .expect("partitioned contig");
            assert_eq!(hash_minimizers, partitioned_minimizers);

            for minimizer in hash_minimizers {
                assert_eq!(
                    hash_sketch.index.get(&minimizer.hash),
                    partitioned_sketch.index.get(&minimizer.hash)
                );
            }
        }

        drop(hash_sketch);
        drop(partitioned_sketch);
        drop(name_sidecar);
        fs::remove_dir_all(directory)?;

        Ok(())
    }

    #[test]
    fn partitioned_sketch_matches_hash_for_odd_and_even_key_counts() -> io::Result<()> {
        let odd_sequence: Vec<u8> = (0..9000)
            .map(|i| b"ACGTGCAATTCG"[i % b"ACGTGCAATTCG".len()])
            .collect();
        assert_partitioned_matches_hash(&odd_sequence, 1, "odd")?;

        let mut state: u64 = 0x9e37_79b9_7f4a_7c15;
        let even_sequence: Vec<u8> = (0..9000)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state & 3) as usize]
            })
            .collect();
        assert_partitioned_matches_hash(&even_sequence, 0, "even")
    }
}
