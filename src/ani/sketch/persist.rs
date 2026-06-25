//! Saving and loading `ReferenceSketch` caches (the `.fasketch` on-disk format).

use std::{
    fs, io,
    io::{BufWriter, Read, Write},
    mem::{align_of, size_of},
    path::Path,
    sync::Arc,
    time::Instant,
};

use crate::ani::{
    align_up, check_memory_limit, checked_section_end, decompress_to_scratch, emit_progress,
    is_gzip_path, load_contig_name_sidecar, slice_as_bytes, write_contig_name_sidecar,
    write_padding, CachedReferenceMetadata, ContigRecord, MinimizerKey, MmapFile,
    MmapReferenceContigs, MmapReferenceIndex, ReferenceContigName, ReferenceContigs, ReferenceFile,
    ReferenceHitMap, ReferenceIndex, ReferenceMinimizer, ReferenceSketch, RuntimeOptions,
    ScratchFile, SeedHit, SketchOutput, SketchParams, SKETCH_KEY_MODE,
    SKETCH_KEY_PACK_PROGRESS_INTERVAL, SKETCH_MAGIC, SKETCH_VERSION,
};
#[cfg(test)]
use crate::ani::{memory_mib, sketch_reference_name};
use boomphf::Mphf;

impl ReferenceSketch {
    /// Save the in-memory reference index as a zero-copy-loadable sketch cache.
    #[cfg(test)]
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn save(
        &self,
        path: &Path,
        kmer_size: usize,
        window_size: usize,
        fragment_length: u32,
        min_fragment_length: u32,
        split_n_run: usize,
        tmp_dir: Option<&Path>,
        runtime_options: RuntimeOptions,
    ) -> io::Result<()> {
        let ReferenceIndex::Hash(index) = &self.index else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "reference sketch cache can only be saved from an in-memory HashMap index",
            ));
        };
        let ReferenceContigs::Owned(contigs) = &self.contigs else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "reference sketch cache can only be saved from owned reference contigs",
            ));
        };

        if let Some(parent) = path
            .parent()
            .filter(|parent| !parent.as_os_str().is_empty())
        {
            fs::create_dir_all(parent)?;
        }

        let save_start: Instant = Instant::now();
        if runtime_options.progress_enabled {
            emit_progress(
                "sketch_save",
                &format!(
                    "event=start\tunique_minimizers={}\tcontigs={}\tfiles={}",
                    index.len(),
                    contigs.len(),
                    self.files.len()
                ),
                save_start,
            );
        }
        check_memory_limit("sketch save start", runtime_options)?;

        let keys: Vec<MinimizerKey> = index.keys().copied().collect::<Vec<_>>();
        if runtime_options.progress_enabled {
            emit_progress(
                "sketch_save",
                &format!("event=keys_collected\tkey_count={}", keys.len()),
                save_start,
            );
        }
        check_memory_limit("sketch save after collecting keys", runtime_options)?;
        let mphf: Mphf<MinimizerKey> = Mphf::new_parallel(1.7, &keys, None);
        if runtime_options.progress_enabled {
            emit_progress(
                "sketch_save",
                &format!("event=mphf_built\tkey_count={}", keys.len()),
                save_start,
            );
        }
        check_memory_limit("sketch save after building MPH", runtime_options)?;
        let mut slot_keys: Vec<MinimizerKey> = vec![0; keys.len()];
        let mut hit_offsets: Vec<u32> = vec![0u32; keys.len()];
        let mut hit_counts: Vec<u32> = vec![0u32; keys.len()];
        let total_hits: usize = index.values().map(Vec::len).sum::<usize>();
        let mut hit_payloads: Vec<SeedHit> = Vec::with_capacity(total_hits);
        let total_reference_minimizers: usize =
            contigs.iter().map(|contig| contig.minimizers.len()).sum();
        let mut contig_records: Vec<ContigRecord> = Vec::with_capacity(contigs.len());
        let (reference_minimizer_scratch, reference_minimizer_file): (ScratchFile, fs::File) =
            ScratchFile::create(tmp_dir, "reference-minimizers")?;
        let mut reference_minimizer_writer: BufWriter<fs::File> =
            BufWriter::new(reference_minimizer_file);
        let mut reference_minimizer_count: usize = 0usize;
        let files: Vec<ReferenceFile> = self
            .files
            .iter()
            .map(|file| ReferenceFile {
                path: sketch_reference_name(&file.path),
                mapped_length: file.mapped_length,
            })
            .collect();

        if runtime_options.progress_enabled {
            #[cfg(debug_assertions)]
            let estimated_pack_bytes: usize = keys
                .len()
                .saturating_mul(size_of::<MinimizerKey>() + size_of::<u64>() + size_of::<u32>())
                .saturating_add(total_hits.saturating_mul(size_of::<SeedHit>()))
                .saturating_add(contigs.len().saturating_mul(size_of::<ContigRecord>()))
                .saturating_add(
                    total_reference_minimizers.saturating_mul(size_of::<ReferenceMinimizer>()),
                );
            #[cfg(debug_assertions)]
            let progress_message: String = format!(
                "event=arrays_allocated\tkey_count={}\thit_count={total_hits}\treference_minimizers={total_reference_minimizers}\testimated_pack_mib={:.3}\ttmp={}",
                keys.len(),
                memory_mib(estimated_pack_bytes),
                reference_minimizer_scratch.path.display()
            );
            #[cfg(not(debug_assertions))]
            let progress_message: String = format!(
                "event=arrays_allocated\tkey_count={}\thit_count={total_hits}\treference_minimizers={total_reference_minimizers}\ttmp={}",
                keys.len(),
                reference_minimizer_scratch.path.display()
            );
            emit_progress("sketch_save", &progress_message, save_start);
        }
        check_memory_limit("sketch save after allocating arrays", runtime_options)?;

        for (key_index, (key, hits)) in index.iter().enumerate() {
            let slot: usize = mphf.hash(key) as usize;
            slot_keys[slot] = *key;
            hit_offsets[slot] = u32::try_from(hit_payloads.len()).map_err(|_| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "shard has more than u32::MAX hit-payload entries ({}); \
                         reduce --shard-minimizers below 4 294 967 295",
                        hit_payloads.len()
                    ),
                )
            })?;
            hit_counts[slot] = hits.len() as u32;
            hit_payloads.extend_from_slice(hits);

            let keys_done: usize = key_index + 1;
            if keys_done.is_multiple_of(SKETCH_KEY_PACK_PROGRESS_INTERVAL)
                || keys_done == index.len()
            {
                if runtime_options.progress_enabled {
                    emit_progress(
                        "sketch_save",
                        &format!(
                            "event=pack_index\tkeys_done={keys_done}\tkey_count={}\thits_done={}",
                            index.len(),
                            hit_payloads.len()
                        ),
                        save_start,
                    );
                }
                check_memory_limit("sketch save while packing index", runtime_options)?;
            }
        }

        for (contig_index, contig) in contigs.iter().enumerate() {
            let minimizer_offset: u64 = reference_minimizer_count as u64;
            let minimizer_count: u32 = u32::try_from(contig.minimizers.len()).map_err(|err| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("reference contig has too many minimizers for sketch cache: {err}"),
                )
            })?;
            let file_id: u32 = u32::try_from(contig.file_id).map_err(|err| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("reference file id exceeds sketch cache limit: {err}"),
                )
            })?;

            contig_records.push(ContigRecord {
                minimizer_offset,
                file_id,
                minimizer_count,
            });
            reference_minimizer_writer.write_all(slice_as_bytes(&contig.minimizers))?;
            reference_minimizer_count += contig.minimizers.len();

            let contigs_done: usize = contig_index + 1;
            if contigs_done.is_multiple_of(crate::ani::SKETCH_CONTIG_PACK_PROGRESS_INTERVAL)
                || contigs_done == contigs.len()
            {
                if runtime_options.progress_enabled {
                    emit_progress(
                        "sketch_save",
                        &format!(
                            "event=pack_contigs\tcontigs_done={contigs_done}\tcontig_count={}\treference_minimizers_done={}",
                            contigs.len(),
                            reference_minimizer_count
                        ),
                        save_start,
                    );
                }
                check_memory_limit("sketch save while packing contigs", runtime_options)?;
            }
        }

        let contig_names: Vec<ReferenceContigName> =
            self.contig_names.clone().unwrap_or_else(|| {
                contigs
                    .iter()
                    .enumerate()
                    .map(|(contig_id, contig)| ReferenceContigName {
                        file_id: contig.file_id,
                        name: format!("contig_{contig_id}"),
                        segment_start: 0,
                        segment_end: 0,
                    })
                    .collect()
            });
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
        write_contig_name_sidecar(
            path,
            &files,
            &contig_names,
            runtime_options.effective_worker_threads(),
        )?;

        let metadata: CachedReferenceMetadata = CachedReferenceMetadata {
            version: SKETCH_VERSION,
            k: kmer_size,
            w: window_size,
            key_mode: SKETCH_KEY_MODE.to_string(),
            fragment_length,
            min_fragment_length,
            split_n_run,
            dust_enabled: false,
            files,
            mphf,
            key_count: slot_keys.len(),
            hit_count: hit_payloads.len(),
            contig_count: contig_records.len(),
            reference_minimizer_count,
        };
        let metadata_bytes: Vec<u8> = serde_json::to_vec(&metadata).map_err(|err| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("failed to encode reference sketch metadata: {err}"),
            )
        })?;
        if runtime_options.progress_enabled {
            emit_progress(
                "sketch_save",
                &format!(
                    "event=metadata_encoded\tmetadata_bytes={}",
                    metadata_bytes.len()
                ),
                save_start,
            );
        }
        check_memory_limit("sketch save after metadata encode", runtime_options)?;
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
            checked_section_end(
                hit_payloads_offset,
                hit_payloads.len(),
                size_of::<SeedHit>(),
            )?,
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

        let mut sketch_output: SketchOutput = SketchOutput::create(
            path,
            tmp_dir,
            false,
            runtime_options.effective_worker_threads(),
        )?;
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
        writer.write_all(slice_as_bytes(&hit_payloads))?;
        write_padding(
            &mut writer,
            contig_records_offset
                - checked_section_end(
                    hit_payloads_offset,
                    hit_payloads.len(),
                    size_of::<SeedHit>(),
                )?,
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
        reference_minimizer_writer.flush()?;
        drop(reference_minimizer_writer);
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
        let mut reference_minimizer_reader: fs::File =
            fs::File::open(&reference_minimizer_scratch.path)?;
        io::copy(&mut reference_minimizer_reader, &mut writer)?;
        let output_file_bytes: u64 = sketch_output.finish()?;
        if runtime_options.progress_enabled {
            emit_progress(
                "sketch_save",
                &format!(
                    "event=complete\tpath={}\tfile_bytes={}",
                    path.display(),
                    output_file_bytes
                ),
                save_start,
            );
        }
        check_memory_limit("sketch save complete", runtime_options)
    }

    /// Load a previously saved sketch cache and mmap its hit arrays.
    pub(crate) fn load(
        path: &Path,
        params: SketchParams,
        load_contig_names: bool,
        tmp_dir: Option<&Path>,
        runtime_options: RuntimeOptions,
    ) -> io::Result<Self> {
        let SketchParams {
            kmer_size,
            window_size,
            fragment_length,
            min_fragment_length,
            split_n_run,
        } = params;
        let load_start: Instant = Instant::now();
        if runtime_options.progress_enabled {
            emit_progress(
                "sketch_load",
                &format!("event=start\tpath={}", path.display()),
                load_start,
            );
        }
        check_memory_limit("sketch load start", runtime_options)?;

        let decompressed_sketch: Option<ScratchFile> = if is_gzip_path(path) {
            if runtime_options.progress_enabled {
                emit_progress(
                    "sketch_load",
                    &format!("event=decompress_start\tpath={}", path.display()),
                    load_start,
                );
            }
            Some(decompress_to_scratch(path, tmp_dir, "decompressed-sketch")?)
        } else {
            None
        };
        let mmap_path: &Path = decompressed_sketch
            .as_ref()
            .map(|scratch| scratch.path.as_path())
            .unwrap_or(path);

        let mut file: fs::File = fs::File::open(mmap_path)?;
        let mut magic: [u8; 8] = [0u8; 8];
        file.read_exact(&mut magic)?;
        if &magic != SKETCH_MAGIC {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "reference sketch cache {} has an invalid magic header",
                    path.display()
                ),
            ));
        }

        let mut metadata_len: [u8; 8] = [0u8; 8];
        file.read_exact(&mut metadata_len)?;
        let metadata_len: usize = u64::from_le_bytes(metadata_len) as usize;
        let metadata_start: usize = SKETCH_MAGIC.len() + size_of::<u64>();
        let metadata_end: usize = metadata_start.checked_add(metadata_len).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "metadata length overflow")
        })?;

        let mmap: Arc<MmapFile> = Arc::new(MmapFile::open(mmap_path)?);
        let bytes: &[u8] = mmap.as_slice();
        if metadata_end > bytes.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "reference sketch cache metadata extends past end of file",
            ));
        }

        let cached: CachedReferenceMetadata =
            serde_json::from_slice(&bytes[metadata_start..metadata_end]).map_err(|err| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("failed to read reference sketch cache: {err}"),
                )
            })?;
        if runtime_options.progress_enabled {
            emit_progress(
                "sketch_load",
                &format!(
                    "event=metadata_loaded\tfiles={}\tcontigs={}\tkey_count={}\thit_count={}\treference_minimizers={}",
                    cached.files.len(),
                    cached.contig_count,
                    cached.key_count,
                    cached.hit_count,
                    cached.reference_minimizer_count
                ),
                load_start,
            );
        }
        check_memory_limit("sketch load after metadata", runtime_options)?;

        if cached.dust_enabled {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "reference sketch cache is incompatible: it was built with the removed --dust filter",
            ));
        }

        if cached.version != SKETCH_VERSION
            || cached.k != kmer_size
            || cached.w != window_size
            || cached.key_mode != SKETCH_KEY_MODE
            || cached.fragment_length != fragment_length
            || cached.min_fragment_length != min_fragment_length
            || cached.split_n_run != split_n_run
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "reference sketch cache is incompatible: version={} k={} w={} key_mode={} fragment_length={} min_fragment_length={} split_n_run={}",
                    cached.version,
                    cached.k,
                    cached.w,
                    cached.key_mode,
                    cached.fragment_length,
                    cached.min_fragment_length,
                    cached.split_n_run
                ),
            ));
        }

        let slot_keys_offset: usize = align_up(metadata_end, 8);
        let hit_offsets_offset: usize = align_up(
            checked_section_end(
                slot_keys_offset,
                cached.key_count,
                size_of::<MinimizerKey>(),
            )?,
            align_of::<u32>(),
        );
        let hit_counts_offset: usize =
            checked_section_end(hit_offsets_offset, cached.key_count, size_of::<u32>())?;
        let hit_payloads_offset: usize = align_up(
            checked_section_end(hit_counts_offset, cached.key_count, size_of::<u32>())?,
            align_of::<SeedHit>(),
        );
        let contig_records_offset: usize = align_up(
            checked_section_end(hit_payloads_offset, cached.hit_count, size_of::<SeedHit>())?,
            align_of::<ContigRecord>(),
        );
        let reference_minimizers_offset: usize = align_up(
            checked_section_end(
                contig_records_offset,
                cached.contig_count,
                size_of::<ContigRecord>(),
            )?,
            align_of::<ReferenceMinimizer>(),
        );
        let file_end: usize = checked_section_end(
            reference_minimizers_offset,
            cached.reference_minimizer_count,
            size_of::<ReferenceMinimizer>(),
        )?;
        if file_end > bytes.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "reference sketch cache array sections extend past end of file",
            ));
        }

        if runtime_options.progress_enabled {
            emit_progress(
                "sketch_load",
                &format!(
                    "event=complete\tmmap_bytes={}\tfiles={}\tcontigs={}\tunique_minimizers={}",
                    bytes.len(),
                    cached.files.len(),
                    cached.contig_count,
                    cached.key_count
                ),
                load_start,
            );
        }
        check_memory_limit("sketch load complete", runtime_options)?;

        let contig_names: Option<Vec<ReferenceContigName>> = if load_contig_names {
            Some(load_contig_name_sidecar(path, cached.contig_count)?)
        } else {
            None
        };

        Ok(Self {
            files: cached.files,
            contigs: ReferenceContigs::Mmap(MmapReferenceContigs {
                mmap: Arc::clone(&mmap),
                contig_count: cached.contig_count,
                reference_minimizer_count: cached.reference_minimizer_count,
                contig_records_offset,
                reference_minimizers_offset,
            }),
            contig_names,
            index: ReferenceIndex::Mphf(MmapReferenceIndex {
                mphf: cached.mphf,
                mmap,
                key_count: cached.key_count,
                hit_count: cached.hit_count,
                slot_keys_offset,
                hit_offsets_offset,
                hit_counts_offset,
                hit_payloads_offset,
            }),
        })
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn save_streamed_cache(
        path: &Path,
        params: SketchParams,
        files: Vec<ReferenceFile>,
        index: &ReferenceHitMap,
        contig_records: Vec<ContigRecord>,
        contig_names: Vec<ReferenceContigName>,
        reference_minimizer_count: usize,
        reference_minimizer_scratch: &ScratchFile,
        tmp_dir: Option<&Path>,
        bgzip: bool,
        runtime_options: RuntimeOptions,
    ) -> io::Result<()> {
        let SketchParams {
            kmer_size,
            window_size,
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
            emit_progress(
                "sketch_save",
                &format!(
                    "event=start\tmode=streaming\tunique_minimizers={}\tcontigs={}\tfiles={}\ttmp={}",
                    index.len(),
                    contig_records.len(),
                    files.len(),
                    reference_minimizer_scratch.path.display()
                ),
                save_start,
            );
        }
        check_memory_limit("streaming sketch save start", runtime_options)?;

        let keys: Vec<MinimizerKey> = index.keys().copied().collect::<Vec<_>>();
        let mphf: Mphf<MinimizerKey> = Mphf::new_parallel(1.7, &keys, None);
        let mut slot_keys: Vec<MinimizerKey> = vec![0; keys.len()];
        let mut hit_offsets: Vec<u32> = vec![0u32; keys.len()];
        let mut hit_counts: Vec<u32> = vec![0u32; keys.len()];
        let total_hits: usize = index.values().map(Vec::len).sum::<usize>();
        let mut hit_payloads: Vec<SeedHit> = Vec::with_capacity(total_hits);

        if runtime_options.progress_enabled {
            emit_progress(
                "sketch_save",
                &format!(
                    "event=arrays_allocated\tmode=streaming\tkey_count={}\thit_count={total_hits}\treference_minimizers={reference_minimizer_count}",
                    keys.len()
                ),
                save_start,
            );
        }
        check_memory_limit(
            "streaming sketch save after allocating arrays",
            runtime_options,
        )?;

        for (key_index, (key, hits)) in index.iter().enumerate() {
            let slot: usize = mphf.hash(key) as usize;
            slot_keys[slot] = *key;
            hit_offsets[slot] = u32::try_from(hit_payloads.len()).map_err(|_| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "shard has more than u32::MAX hit-payload entries ({}); \
                         reduce --shard-minimizers below 4 294 967 295",
                        hit_payloads.len()
                    ),
                )
            })?;
            hit_counts[slot] = hits.len() as u32;
            hit_payloads.extend_from_slice(hits);

            let keys_done: usize = key_index + 1;
            if keys_done.is_multiple_of(SKETCH_KEY_PACK_PROGRESS_INTERVAL)
                || keys_done == index.len()
            {
                if runtime_options.progress_enabled {
                    emit_progress(
                        "sketch_save",
                        &format!(
                            "event=pack_index\tmode=streaming\tkeys_done={keys_done}\tkey_count={}\thits_done={}",
                            index.len(),
                            hit_payloads.len()
                        ),
                        save_start,
                    );
                }
                check_memory_limit("streaming sketch save while packing index", runtime_options)?;
            }
        }

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
        write_contig_name_sidecar(
            path,
            &files,
            &contig_names,
            runtime_options.effective_worker_threads(),
        )?;

        let metadata: CachedReferenceMetadata = CachedReferenceMetadata {
            version: SKETCH_VERSION,
            k: kmer_size,
            w: window_size,
            key_mode: SKETCH_KEY_MODE.to_string(),
            fragment_length,
            min_fragment_length,
            split_n_run,
            dust_enabled: false,
            files,
            mphf,
            key_count: slot_keys.len(),
            hit_count: hit_payloads.len(),
            contig_count: contig_records.len(),
            reference_minimizer_count,
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
            align_of::<u64>(),
        );
        let hit_counts_offset: usize =
            checked_section_end(hit_offsets_offset, hit_offsets.len(), size_of::<u64>())?;
        let hit_payloads_offset: usize = align_up(
            checked_section_end(hit_counts_offset, hit_counts.len(), size_of::<u32>())?,
            align_of::<SeedHit>(),
        );
        let contig_records_offset: usize = align_up(
            checked_section_end(
                hit_payloads_offset,
                hit_payloads.len(),
                size_of::<SeedHit>(),
            )?,
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

        let mut sketch_output: SketchOutput = SketchOutput::create(
            path,
            tmp_dir,
            bgzip,
            runtime_options.effective_worker_threads(),
        )?;
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
        writer.write_all(slice_as_bytes(&hit_payloads))?;
        write_padding(
            &mut writer,
            contig_records_offset
                - checked_section_end(
                    hit_payloads_offset,
                    hit_payloads.len(),
                    size_of::<SeedHit>(),
                )?,
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
            emit_progress(
                "sketch_save",
                &format!(
                    "event=complete\tmode=streaming\tpath={}\tfile_bytes={}",
                    path.display(),
                    output_file_bytes
                ),
                save_start,
            );
        }
        check_memory_limit("streaming sketch save complete", runtime_options)
    }
}

#[cfg(test)]
mod tests {
    use crate::ani::{
        contig_sidecar_path, ReferenceContig, ReferenceContigName, ReferenceContigs, ReferenceFile,
        ReferenceHitMap, ReferenceIndex, ReferenceMinimizer, ReferenceSketch, RuntimeOptions,
        SeedHit, SketchParams, DEFAULT_FRAGMENT_LENGTH, DEFAULT_KMER_SIZE,
        DEFAULT_MIN_FRAGMENT_LENGTH, DEFAULT_WINDOW_SIZE,
    };
    use std::{env, fs, io, path::PathBuf, time::Instant};

    #[test]
    fn loaded_sketch_contig_minimizers_match_owned_minimizers() -> io::Result<()> {
        let contigs: Vec<ReferenceContig> = vec![
            ReferenceContig {
                file_id: 0,
                minimizers: vec![
                    ReferenceMinimizer {
                        hash: 11,
                        position: 3,
                    },
                    ReferenceMinimizer {
                        hash: 17,
                        position: 9,
                    },
                ],
            },
            ReferenceContig {
                file_id: 0,
                minimizers: vec![ReferenceMinimizer {
                    hash: 23,
                    position: 4,
                }],
            },
        ];
        let mut index: ReferenceHitMap = ReferenceHitMap::default();
        for (contig_id, contig) in contigs.iter().enumerate() {
            for minimizer in &contig.minimizers {
                index.entry(minimizer.hash).or_default().push(SeedHit {
                    reference_contig_id: contig_id as u32,
                    position: minimizer.position,
                });
            }
        }

        let sketch: ReferenceSketch = ReferenceSketch {
            files: vec![ReferenceFile {
                path: "/tmp/ref.fa".to_string(),
                mapped_length: 3000,
            }],
            contigs: ReferenceContigs::Owned(contigs.clone()),
            contig_names: Some(vec![
                ReferenceContigName {
                    file_id: 0,
                    name: "ref_contig_a".to_string(),
                    segment_start: 0,
                    segment_end: 100,
                },
                ReferenceContigName {
                    file_id: 0,
                    name: "ref_contig_b".to_string(),
                    segment_start: 200,
                    segment_end: 300,
                },
            ]),
            index: ReferenceIndex::Hash(index),
        };
        let path: PathBuf = env::temp_dir().join(format!(
            "fasterani_mmap_contigs_{}_{}.fasketch",
            std::process::id(),
            Instant::now().elapsed().as_nanos()
        ));

        sketch.save(
            &path,
            DEFAULT_KMER_SIZE,
            DEFAULT_WINDOW_SIZE,
            DEFAULT_FRAGMENT_LENGTH,
            DEFAULT_MIN_FRAGMENT_LENGTH,
            0,
            None,
            RuntimeOptions::default(),
        )?;
        let loaded: ReferenceSketch = ReferenceSketch::load(
            &path,
            SketchParams {
                kmer_size: DEFAULT_KMER_SIZE,
                window_size: DEFAULT_WINDOW_SIZE,
                fragment_length: DEFAULT_FRAGMENT_LENGTH,
                min_fragment_length: DEFAULT_MIN_FRAGMENT_LENGTH,
                split_n_run: 0,
            },
            true,
            None,
            RuntimeOptions::default(),
        )?;
        let loaded_contig_names: &[ReferenceContigName] =
            loaded.contig_names.as_deref().expect("loaded sidecar");
        assert_eq!(loaded_contig_names[0].name, "ref_contig_a");
        assert_eq!(loaded_contig_names[1].segment_start, 200);
        fs::remove_file(&path)?;
        fs::remove_file(contig_sidecar_path(&path))?;

        assert_eq!(loaded.files[0].path, "ref.fa");
        assert_eq!(loaded.contigs.len(), contigs.len());
        for (contig_id, contig) in contigs.iter().enumerate() {
            assert_eq!(loaded.contigs.file_id(contig_id), Some(contig.file_id));
            assert_eq!(
                loaded.contigs.minimizers(contig_id).expect("loaded contig"),
                contig.minimizers.as_slice()
            );
        }

        Ok(())
    }
}
