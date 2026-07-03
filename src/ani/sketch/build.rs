//! Building an in-memory `ReferenceSketch` from reference FASTA files.

#[cfg(debug_assertions)]
use std::mem::size_of;
use std::{io, time::Instant};

use crate::ani::{
    canonical_minimizers_with_positions, emit_progress, mapped_length_from_fragment_ranges,
    open_fasta_reader, split_sequence_ranges, FastaInput, ReferenceContig, ReferenceContigName,
    ReferenceContigs, ReferenceFile, ReferenceHitMap, ReferenceIndex, ReferenceMinimizer,
    ReferenceSketch, RuntimeOptions, SeedHit, SketchParams, REFERENCE_PROGRESS_INTERVAL,
};
#[cfg(debug_assertions)]
use crate::ani::{
    memory_mib, reference_build_struct_bytes, ContigRecord, MinimizerKey, ReferenceMemoryEstimate,
};
use noodles::fasta;

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
        let SketchParams {
            kmer_size,
            window_size,
            fragment_length,
            min_fragment_length,
            split_n_run,
        } = params;
        let mut files: Vec<ReferenceFile> = Vec::new();
        let mut contigs: Vec<ReferenceContig> = Vec::new();
        let mut index: ReferenceHitMap = ReferenceHitMap::default();
        let build_start: Instant = Instant::now();
        let mut total_reference_minimizers: usize = 0usize;
        let mut total_seed_hits: usize = 0usize;
        let mut contig_names: Vec<ReferenceContigName> = Vec::new();
        files.reserve(references.len());
        contigs.reserve(references.len());

        if runtime_options.progress_enabled {
            emit_progress(
                "reference_build",
                &format!(
                    "event=start\tfiles_total={}\tsplit_n_run={split_n_run}",
                    references.len()
                ),
                build_start,
            );
        }
        for (file_id, reference) in references.iter().enumerate() {
            let mut reader: fasta::io::Reader<Box<dyn io::BufRead>> =
                open_fasta_reader(&reference.open)?;
            let mut mapped_length: u64 = 0u64;

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

                    let reference_contig_id: usize = contigs.len();
                    let sketch_capacity: usize =
                        segment_sequence.len() / (window_size / 2).max(1) + 1;
                    let mut reference_minimizers: Vec<ReferenceMinimizer> =
                        Vec::with_capacity(sketch_capacity);

                    if segment_sequence.len() >= kmer_size && segment_sequence.len() >= window_size
                    {
                        for (hash, position) in canonical_minimizers_with_positions(
                            segment_sequence,
                            kmer_size,
                            window_size,
                        ) {
                            reference_minimizers.push(ReferenceMinimizer { hash, position });
                        }
                    }

                    reference_minimizers.sort_unstable_by_key(|minimizer| minimizer.position);
                    index.reserve(reference_minimizers.len());

                    for minimizer in &reference_minimizers {
                        index.entry(minimizer.hash).or_default().push(SeedHit {
                            reference_contig_id: reference_contig_id as u32,
                            position: minimizer.position,
                        });
                    }

                    total_reference_minimizers += reference_minimizers.len();
                    total_seed_hits += reference_minimizers.len();

                    contigs.push(ReferenceContig {
                        file_id,
                        minimizers: reference_minimizers,
                    });
                    contig_names.push(ReferenceContigName {
                        file_id,
                        name: record_name.clone(),
                        segment_start,
                        segment_end,
                    });
                }
            }

            files.push(ReferenceFile {
                path: reference.label.clone(),
                mapped_length,
            });

            let files_done: usize = file_id + 1;
            if runtime_options.progress_enabled
                && (files_done.is_multiple_of(REFERENCE_PROGRESS_INTERVAL)
                    || files_done == references.len())
            {
                #[cfg(debug_assertions)]
                let estimated_struct_bytes: usize = reference_build_struct_bytes(
                    total_reference_minimizers,
                    total_seed_hits,
                    index.len(),
                );
                #[cfg(debug_assertions)]
                let progress_message: String = format!(
                    "event=files\tfiles_done={files_done}\tfiles_total={}\tcontigs={}\treference_minimizers={total_reference_minimizers}\tunique_minimizers={}\tseed_hits={total_seed_hits}\testimated_struct_mib={:.3}",
                    references.len(),
                    contigs.len(),
                    index.len(),
                    memory_mib(estimated_struct_bytes)
                );
                #[cfg(not(debug_assertions))]
                let progress_message: String = format!(
                    "event=files\tfiles_done={files_done}\tfiles_total={}\tcontigs={}\treference_minimizers={total_reference_minimizers}\tunique_minimizers={}\tseed_hits={total_seed_hits}",
                    references.len(),
                    contigs.len(),
                    index.len()
                );
                emit_progress("reference_build", &progress_message, build_start);
            }
        }

        if runtime_options.progress_enabled {
            emit_progress(
                "reference_build",
                &format!(
                    "event=complete\tfiles_done={}\tcontigs={}\treference_minimizers={total_reference_minimizers}\tunique_minimizers={}\tseed_hits={total_seed_hits}",
                    references.len(),
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
