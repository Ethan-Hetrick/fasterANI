//! Building an in-memory `ReferenceSketch` from reference FASTA files.

#[cfg(debug_assertions)]
use std::mem::size_of;
use std::{
    io,
    sync::atomic::{AtomicUsize, Ordering as AtomicOrdering},
    time::Instant,
};

use crate::ani::{
    emit_progress, for_each_extracted_reference_segment, FastaInput, MinimizerKey, ReferenceContig,
    ReferenceContigName, ReferenceContigs, ReferenceFile, ReferenceIndex, ReferenceMinimizer,
    ReferenceSketch, RuntimeOptions, SeedHit, SketchParams, TransientReferenceIndex,
    REFERENCE_PROGRESS_INTERVAL,
};
#[cfg(debug_assertions)]
use crate::ani::{memory_mib, reference_build_struct_bytes, ContigRecord, ReferenceMemoryEstimate};
use rayon::prelude::*;

struct FileBuild {
    file: ReferenceFile,
    contigs: Vec<ReferenceContig>,
    contig_names: Vec<ReferenceContigName>,
    keyed_hits: Vec<(MinimizerKey, SeedHit)>,
    minimizer_count: usize,
}

fn collect_reference_file(reference: &FastaInput, params: SketchParams) -> io::Result<FileBuild> {
    let mut contigs: Vec<ReferenceContig> = Vec::new();
    let mut contig_names: Vec<ReferenceContigName> = Vec::new();
    let mut keyed_hits: Vec<(MinimizerKey, SeedHit)> = Vec::new();
    let mut minimizer_count: usize = 0usize;

    let mapped_length: u64 = for_each_extracted_reference_segment(reference, params, |segment| {
        let local_contig_id_u32: u32 = u32::try_from(contigs.len()).map_err(|err| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "reference {} record {:?}: contig id exceeds u32: {err}",
                    reference.label, segment.record_name
                ),
            )
        })?;
        keyed_hits.reserve(segment.minimizers.len());

        for minimizer in &segment.minimizers {
            keyed_hits.push((
                minimizer.hash,
                SeedHit {
                    reference_contig_id: local_contig_id_u32,
                    position: minimizer.position,
                },
            ));
        }

        minimizer_count = minimizer_count
            .checked_add(segment.minimizers.len())
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "reference {} minimizer count exceeds usize",
                        reference.label
                    ),
                )
            })?;
        contigs.push(ReferenceContig {
            file_id: 0,
            minimizers: segment.minimizers,
        });
        contig_names.push(ReferenceContigName {
            file_id: 0,
            name: segment.record_name,
            segment_start: segment.segment_start,
            segment_end: segment.segment_end,
        });
        Ok(())
    })?;

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

fn unique_key_count(hits: &[(MinimizerKey, SeedHit)]) -> usize {
    if hits.is_empty() {
        0
    } else {
        1 + hits
            .windows(2)
            .filter(|pair| pair[0].0 != pair[1].0)
            .count()
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
            ReferenceIndex::HashShards(shards) => {
                let seed_hits: usize = shards
                    .iter()
                    .flat_map(|shard| shard.index.values())
                    .map(Vec::len)
                    .sum();
                let seed_hit_vec_bytes: usize = shards
                    .iter()
                    .flat_map(|shard| shard.index.values())
                    .map(|hits| hits.capacity() * size_of::<SeedHit>())
                    .sum();
                let unique_index_keys: usize = shards.iter().map(|shard| shard.index.len()).sum();
                let hash_index_rough_bytes: usize = shards
                    .iter()
                    .map(|shard| {
                        shard.index.capacity()
                            * (size_of::<MinimizerKey>() + size_of::<Vec<SeedHit>>() + 8)
                    })
                    .sum();

                ReferenceMemoryEstimate {
                    reference_minimizers,
                    reference_minimizer_vec_bytes,
                    unique_index_keys,
                    seed_hits,
                    seed_hit_vec_bytes,
                    hash_index_rough_bytes,
                    ..ReferenceMemoryEstimate::default()
                }
            }
            ReferenceIndex::Transient(index) => ReferenceMemoryEstimate {
                reference_minimizers,
                reference_minimizer_vec_bytes,
                unique_index_keys: index.keys.len(),
                seed_hits: index.hit_payloads.len(),
                seed_hit_vec_bytes: index.hit_payloads.capacity() * size_of::<SeedHit>(),
                hash_index_rough_bytes: index.keys.capacity() * size_of::<MinimizerKey>()
                    + index.hit_offsets.capacity() * size_of::<usize>(),
                ..ReferenceMemoryEstimate::default()
            },
            ReferenceIndex::Mphf(index) => ReferenceMemoryEstimate {
                reference_minimizers,
                reference_minimizer_vec_bytes,
                unique_index_keys: index.key_count,
                seed_hits: index.hit_count,
                mmap_file_bytes: index.mmap.len,
                mmap_slot_key_bytes: index.key_count * size_of::<MinimizerKey>(),
                mmap_hit_offset_bytes: index.key_count * size_of::<u32>(),
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
                &format!("event=start\tfiles_total={total_files}\tsplit_n_run={split_n_run}"),
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

        let flatten_start: Instant = Instant::now();
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

        if let Some(start) = merge_start {
            emit_progress(
                "reference_merge",
                &format!(
                    "event=flatten_complete\tphase_ms={:.3}\thit_records={}\tcontigs={}\tfiles={}\tworker_threads={worker_threads}",
                    flatten_start.elapsed().as_secs_f64() * 1000.0,
                    all_hits.len(),
                    contigs.len(),
                    files.len()
                ),
                start,
            );
        }

        let sort_start: Instant = Instant::now();
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
        let sorted_unique_minimizers: usize = unique_key_count(&all_hits);

        if let Some(start) = merge_start {
            emit_progress(
                "reference_merge",
                &format!(
                    "event=sort_complete\tphase_ms={:.3}\thit_records={}\tunique_minimizers={}\tworker_threads={worker_threads}",
                    sort_start.elapsed().as_secs_f64() * 1000.0,
                    all_hits.len(),
                    sorted_unique_minimizers
                ),
                start,
            );
        }

        let index_build_start: Instant = Instant::now();
        let index: ReferenceIndex =
            ReferenceIndex::Transient(TransientReferenceIndex::from_sorted_hits_with_key_count(
                &all_hits,
                sorted_unique_minimizers,
            ));
        let target_ranges: usize = usize::from(!all_hits.is_empty());
        let actual_ranges: usize = target_ranges;
        let min_range_hits: usize = if all_hits.is_empty() {
            0
        } else {
            all_hits.len()
        };
        let max_range_hits: usize = min_range_hits;
        let index_mode: &'static str = "transient_flat";
        let shards: usize = 0;

        if let Some(start) = merge_start {
            emit_progress(
                "reference_merge",
                &format!(
                    "event=index_build_complete\tphase_ms={:.3}\thit_records={}\tunique_minimizers={}\tworker_threads={worker_threads}\tindex_mode={index_mode}\tshards={shards}",
                    index_build_start.elapsed().as_secs_f64() * 1000.0,
                    all_hits.len(),
                    index.len(),
                ),
                start,
            );
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
                    "event=complete\thit_records={}\tcontigs={}\tunique_minimizers={}\tseed_hits={total_seed_hits}\tworker_threads={worker_threads}\ttarget_ranges={target_ranges}\tactual_ranges={actual_ranges}\tmin_range_hits={min_range_hits}\tmax_range_hits={max_range_hits}\tindex_mode={index_mode}\tshards={shards}\testimated_struct_mib={:.3}",
                    all_hits.len(),
                    contigs.len(),
                    index.len(),
                    memory_mib(estimated_struct_bytes)
                )
            };
            #[cfg(not(debug_assertions))]
            let progress_message: String = format!(
                "event=complete\thit_records={}\tcontigs={}\tunique_minimizers={}\tseed_hits={total_seed_hits}\tworker_threads={worker_threads}\ttarget_ranges={target_ranges}\tactual_ranges={actual_ranges}\tmin_range_hits={min_range_hits}\tmax_range_hits={max_range_hits}\tindex_mode={index_mode}\tshards={shards}",
                all_hits.len(),
                contigs.len(),
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
            index,
            global_frequencies: None,
        })
    }
}

#[cfg(test)]
mod tests {
    use crate::ani::{
        FastaInput, ReferenceIndex, ReferenceSketch, RuntimeOptions, SketchParams,
        DEFAULT_FRAGMENT_LENGTH, DEFAULT_KMER_SIZE, DEFAULT_MINIMIZER_HASH_SEED,
        DEFAULT_MIN_FRAGMENT_LENGTH, DEFAULT_SPLIT_N_RUN, DEFAULT_WINDOW_SIZE,
    };
    use std::{env, fs, io, path::PathBuf, time::Instant};

    #[test]
    fn no_save_collect_builds_transient_reference_index() -> io::Result<()> {
        let sequence: String = (0..300)
            .map(|i| b"ACGTGCAATTCG"[i % b"ACGTGCAATTCG".len()] as char)
            .collect();
        let reference_path: PathBuf = env::temp_dir().join(format!(
            "fasterani_transient_ref_{}_{}.fa",
            std::process::id(),
            Instant::now().elapsed().as_nanos()
        ));
        fs::write(&reference_path, format!(">ref\n{sequence}\n"))?;
        let references: Vec<FastaInput> = vec![FastaInput::from_path(
            reference_path.to_string_lossy().into_owned(),
        )];

        let sketch: ReferenceSketch = ReferenceSketch::collect(
            &references,
            SketchParams {
                kmer_size: DEFAULT_KMER_SIZE,
                window_size: DEFAULT_WINDOW_SIZE,
                minimizer_hash_seed: DEFAULT_MINIMIZER_HASH_SEED,
                fragment_length: DEFAULT_FRAGMENT_LENGTH,
                min_fragment_length: DEFAULT_MIN_FRAGMENT_LENGTH,
                split_n_run: DEFAULT_SPLIT_N_RUN,
            },
            RuntimeOptions::default(),
        )?;

        let first_minimizer = sketch
            .contigs
            .minimizers(0)
            .and_then(|minimizers| minimizers.first())
            .copied()
            .expect("reference minimizer");
        match &sketch.index {
            ReferenceIndex::Transient(index) => {
                let hits = index
                    .get(&first_minimizer.hash)
                    .expect("transient minimizer hit");
                assert!(hits.iter().any(|hit| {
                    hit.reference_contig_id == 0 && hit.position == first_minimizer.position
                }));
            }
            _ => panic!("no-save collect should build transient index"),
        }

        fs::remove_file(reference_path)?;
        Ok(())
    }
}
