//! Building an in-memory `ReferenceSketch` from reference FASTA files.

#[cfg(debug_assertions)]
use std::mem::size_of;
use std::{
    io,
    sync::atomic::{AtomicUsize, Ordering as AtomicOrdering},
    time::Instant,
};

use crate::ani::{
    canonical_minimizers_with_positions, emit_progress, mapped_length_from_fragment_ranges,
    open_fasta_reader, split_sequence_ranges, FastaInput, MinimizerKey, ReferenceContig,
    ReferenceContigName, ReferenceContigs, ReferenceFile, ReferenceHitMap, ReferenceIndex,
    ReferenceMinimizer, ReferenceSketch, RuntimeOptions, SeedHit, SketchParams,
    REFERENCE_PROGRESS_INTERVAL,
};
#[cfg(debug_assertions)]
use crate::ani::{memory_mib, reference_build_struct_bytes, ContigRecord, ReferenceMemoryEstimate};
use noodles::fasta;
use rayon::prelude::*;

struct FileBuild {
    file: ReferenceFile,
    contigs: Vec<ReferenceContig>,
    contig_names: Vec<ReferenceContigName>,
    keyed_hits: Vec<(MinimizerKey, SeedHit)>,
    minimizer_count: usize,
}

fn collect_reference_file(reference: &FastaInput, params: SketchParams) -> io::Result<FileBuild> {
    let SketchParams {
        kmer_size,
        window_size,
        fragment_length,
        min_fragment_length,
        split_n_run,
    } = params;
    let mut reader: fasta::io::Reader<Box<dyn io::BufRead>> = open_fasta_reader(&reference.open)?;
    let mut mapped_length: u64 = 0u64;
    let mut contigs: Vec<ReferenceContig> = Vec::new();
    let mut contig_names: Vec<ReferenceContigName> = Vec::new();
    let mut keyed_hits: Vec<(MinimizerKey, SeedHit)> = Vec::new();
    let mut minimizer_count: usize = 0usize;

    for result in reader.records() {
        let record: fasta::Record = result.map_err(|err| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "failed to read FASTA record from reference {}: {err}",
                    reference.label
                ),
            )
        })?;
        let record_name: String = String::from_utf8_lossy(record.name()).into_owned();
        let sequence: &fasta::record::Sequence = record.sequence();
        let sequence_bytes: &[u8] = sequence.as_ref();

        for segment_range in split_sequence_ranges(sequence_bytes, split_n_run) {
            let segment_start: u32 = u32::try_from(segment_range.start).map_err(|err| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("reference segment start exceeds u32: {err}"),
                )
            })?;
            let segment_end: u32 = u32::try_from(segment_range.end).map_err(|err| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("reference segment end exceeds u32: {err}"),
                )
            })?;
            let segment_sequence: &[u8] = &sequence_bytes[segment_range];

            mapped_length += mapped_length_from_fragment_ranges(
                segment_sequence.len(),
                fragment_length,
                min_fragment_length,
            );

            let local_contig_id: usize = contigs.len();
            let local_contig_id_u32: u32 = u32::try_from(local_contig_id).map_err(|err| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("reference contig id exceeds u32: {err}"),
                )
            })?;
            let sketch_capacity: usize = segment_sequence.len() / (window_size / 2).max(1) + 1;
            let mut reference_minimizers: Vec<ReferenceMinimizer> =
                Vec::with_capacity(sketch_capacity);

            if segment_sequence.len() >= kmer_size && segment_sequence.len() >= window_size {
                for (hash, position) in
                    canonical_minimizers_with_positions(segment_sequence, kmer_size, window_size)
                {
                    reference_minimizers.push(ReferenceMinimizer { hash, position });
                }
            }

            reference_minimizers.sort_unstable_by_key(|minimizer| minimizer.position);
            keyed_hits.reserve(reference_minimizers.len());

            for minimizer in &reference_minimizers {
                keyed_hits.push((
                    minimizer.hash,
                    SeedHit {
                        reference_contig_id: local_contig_id_u32,
                        position: minimizer.position,
                    },
                ));
            }

            minimizer_count += reference_minimizers.len();
            contigs.push(ReferenceContig {
                file_id: 0,
                minimizers: reference_minimizers,
            });
            contig_names.push(ReferenceContigName {
                file_id: 0,
                name: record_name.clone(),
                segment_start,
                segment_end,
            });
        }
    }

    Ok(FileBuild {
        file: ReferenceFile {
            path: reference.label.clone(),
            mapped_length,
        },
        contigs,
        contig_names,
        keyed_hits,
        minimizer_count,
    })
}

fn emit_reference_files_progress(files_done: usize, total_files: usize, build_start: Instant) {
    if files_done.is_multiple_of(REFERENCE_PROGRESS_INTERVAL) || files_done == total_files {
        emit_progress(
            "reference_build",
            &format!("event=files\tfiles_done={files_done}\tfiles_total={total_files}"),
            build_start,
        );
    }
}

impl ReferenceSketch {
    #[cfg(debug_assertions)]
    pub(crate) fn memory_estimate(&self) -> ReferenceMemoryEstimate {
        let reference_minimizers: usize = self.contigs.total_minimizers();
        let reference_minimizer_vec_bytes: usize = self.contigs.owned_minimizer_capacity_bytes();

        match &self.index {
            ReferenceIndex::Hash(index) => {
                let seed_hits: usize = index.values().map(Vec::len).sum();
                let seed_hit_vec_bytes: usize = index
                    .values()
                    .map(|hits| hits.capacity() * size_of::<SeedHit>())
                    .sum();
                let hash_index_rough_bytes: usize =
                    index.capacity() * (size_of::<MinimizerKey>() + size_of::<Vec<SeedHit>>() + 8);

                ReferenceMemoryEstimate {
                    reference_minimizers,
                    reference_minimizer_vec_bytes,
                    unique_index_keys: index.len(),
                    seed_hits,
                    seed_hit_vec_bytes,
                    hash_index_rough_bytes,
                    ..ReferenceMemoryEstimate::default()
                }
            }
            ReferenceIndex::Mphf(index) => ReferenceMemoryEstimate {
                reference_minimizers,
                reference_minimizer_vec_bytes,
                unique_index_keys: index.key_count,
                seed_hits: index.hit_count,
                mmap_file_bytes: index.mmap.len,
                mmap_slot_key_bytes: index.key_count * size_of::<MinimizerKey>(),
                mmap_hit_offset_bytes: index.key_count * size_of::<u64>(),
                mmap_hit_count_bytes: index.key_count * size_of::<u32>(),
                mmap_hit_payload_bytes: index.hit_count * size_of::<SeedHit>(),
                mmap_contig_record_bytes: self.contigs.len() * size_of::<ContigRecord>(),
                mmap_reference_minimizer_bytes: reference_minimizers
                    * size_of::<ReferenceMinimizer>(),
                ..ReferenceMemoryEstimate::default()
            },
        }
    }

    /// Build a reference sketch from all provided reference FASTA files.
    pub(crate) fn collect(
        references: &[FastaInput],
        params: SketchParams,
        runtime_options: RuntimeOptions,
    ) -> io::Result<Self> {
        let SketchParams { split_n_run, .. } = params;
        let build_start: Instant = Instant::now();
        let total_files: usize = references.len();
        let worker_threads: usize = runtime_options.effective_worker_threads();
        let pool: Option<rayon::ThreadPool> = if worker_threads > 1 {
            Some(
                rayon::ThreadPoolBuilder::new()
                    .num_threads(worker_threads)
                    .build()
                    .map_err(|err| {
                        io::Error::new(
                            io::ErrorKind::InvalidInput,
                            format!("failed to initialize reference build thread pool: {err}"),
                        )
                    })?,
            )
        } else {
            None
        };

        if runtime_options.progress_enabled {
            emit_progress(
                "reference_build",
                &format!(
                    "event=start\tfiles_total={}\tsplit_n_run={split_n_run}",
                    total_files
                ),
                build_start,
            );
        }

        let mut per_file: Vec<FileBuild> = if let Some(pool) = pool.as_ref() {
            let completed_files: AtomicUsize = AtomicUsize::new(0);
            pool.install(|| {
                references
                    .par_iter()
                    .map(|reference| {
                        let build: FileBuild = collect_reference_file(reference, params)?;
                        let files_done: usize =
                            completed_files.fetch_add(1, AtomicOrdering::Relaxed) + 1;
                        if runtime_options.progress_enabled {
                            emit_reference_files_progress(files_done, total_files, build_start);
                        }
                        Ok(build)
                    })
                    .collect::<io::Result<Vec<_>>>()
            })?
        } else {
            let mut builds: Vec<FileBuild> = Vec::with_capacity(total_files);
            for reference in references {
                builds.push(collect_reference_file(reference, params)?);
                if runtime_options.progress_enabled {
                    emit_reference_files_progress(builds.len(), total_files, build_start);
                }
            }
            builds
        };

        let merge_start: Option<Instant> = runtime_options.progress_enabled.then(Instant::now);
        if let Some(start) = merge_start {
            emit_progress("reference_merge", "event=start", start);
        }

        let total_contigs: usize = per_file.iter().map(|build| build.contigs.len()).sum();
        let total_reference_minimizers: usize =
            per_file.iter().map(|build| build.minimizer_count).sum();
        let total_seed_hits: usize = total_reference_minimizers;
        let mut files: Vec<ReferenceFile> = Vec::with_capacity(total_files);
        let mut contigs: Vec<ReferenceContig> = Vec::with_capacity(total_contigs);
        let mut contig_names: Vec<ReferenceContigName> = Vec::with_capacity(total_contigs);
        let mut all_hits: Vec<(MinimizerKey, SeedHit)> = Vec::with_capacity(total_seed_hits);

        for (file_id, mut built) in per_file.drain(..).enumerate() {
            let contig_offset: usize = contigs.len();
            let contig_offset_u32: u32 = u32::try_from(contig_offset).map_err(|err| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("reference contig offset exceeds u32: {err}"),
                )
            })?;

            for contig in &mut built.contigs {
                contig.file_id = file_id;
            }
            for contig_name in &mut built.contig_names {
                contig_name.file_id = file_id;
            }
            for (_key, hit) in &mut built.keyed_hits {
                hit.reference_contig_id = hit
                    .reference_contig_id
                    .checked_add(contig_offset_u32)
                    .ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            "reference contig id exceeds u32",
                        )
                    })?;
            }

            files.push(built.file);
            contigs.extend(built.contigs);
            contig_names.extend(built.contig_names);
            all_hits.extend(built.keyed_hits);
        }

        if let Some(pool) = pool.as_ref() {
            pool.install(|| {
                all_hits.par_sort_unstable_by_key(|&(key, hit)| {
                    (key, hit.reference_contig_id, hit.position)
                });
            });
        } else {
            all_hits
                .sort_unstable_by_key(|&(key, hit)| (key, hit.reference_contig_id, hit.position));
        }

        let unique_minimizer_count: usize = if all_hits.is_empty() {
            0
        } else {
            1 + all_hits
                .windows(2)
                .filter(|pair| pair[0].0 != pair[1].0)
                .count()
        };
        let mut index: ReferenceHitMap = ReferenceHitMap::default();
        index.reserve(unique_minimizer_count);
        let mut cursor: usize = 0usize;
        while cursor < all_hits.len() {
            let key: MinimizerKey = all_hits[cursor].0;
            let start: usize = cursor;
            cursor += 1;

            while cursor < all_hits.len() && all_hits[cursor].0 == key {
                cursor += 1;
            }

            let mut hits: Vec<SeedHit> = Vec::with_capacity(cursor - start);
            hits.extend(all_hits[start..cursor].iter().map(|&(_key, hit)| hit));
            index.insert(key, hits);
        }

        if let Some(start) = merge_start {
            #[cfg(debug_assertions)]
            let progress_message: String = {
                let estimated_struct_bytes: usize = reference_build_struct_bytes(
                    total_reference_minimizers,
                    total_seed_hits,
                    index.len(),
                );
                format!(
                    "event=complete\tunique_minimizers={}\tseed_hits={total_seed_hits}\testimated_struct_mib={:.3}",
                    index.len(),
                    memory_mib(estimated_struct_bytes)
                )
            };
            #[cfg(not(debug_assertions))]
            let progress_message: String = format!(
                "event=complete\tunique_minimizers={}\tseed_hits={total_seed_hits}",
                index.len()
            );
            emit_progress("reference_merge", &progress_message, start);
        }
        if runtime_options.progress_enabled {
            emit_progress(
                "reference_build",
                &format!(
                    "event=complete\tfiles_done={}\tcontigs={}\treference_minimizers={total_reference_minimizers}\tunique_minimizers={}\tseed_hits={total_seed_hits}",
                    total_files,
                    contigs.len(),
                    index.len()
                ),
                build_start,
            );
        }

        Ok(Self {
            files,
            contigs: ReferenceContigs::Owned(contigs),
            contig_names: Some(contig_names),
            index: ReferenceIndex::Hash(index),
        })
    }
}
