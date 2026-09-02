//! Ordered, bounded extraction of reference minimizers from FASTA inputs.

use std::{io, ops::Range};

use noodles::fasta;
use rayon::prelude::*;

use crate::ani::{
    canonical_minimizers_with_positions, canonical_minimizers_with_super_kmers,
    mapped_length_from_fragment_ranges, open_fasta_reader, split_sequence_ranges, FastaInput,
    ReferenceMinimizer, SketchParams,
};

// Keep enough work in flight to make a multi-record FASTA useful to Rayon without retaining an
// entire large input. A single record can exceed the byte threshold, but no second record is then
// admitted to that batch.
const EXTRACT_BATCH_RECORDS: usize = 64;
const EXTRACT_BATCH_BASES: usize = 32 * 1024 * 1024;
const EXTRACT_SEQUENCE_CHUNK_WINDOWS: usize = 2 * 1024 * 1024;
const EXTRACT_CLEAN_RANGE_BATCH: usize = 1024;

pub(crate) struct ExtractedReferenceSegment {
    pub(crate) record_name: String,
    pub(crate) segment_start: u32,
    pub(crate) segment_end: u32,
    pub(crate) minimizers: Vec<ReferenceMinimizer>,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(crate) struct ReferenceExtractionStats {
    pub(crate) mapped_length: u64,
    pub(crate) original_length: u64,
}

struct PendingRecord {
    name: String,
    sequence: Vec<u8>,
}

struct ExtractedRecord {
    mapped_length: u64,
    segments: Vec<ExtractedReferenceSegment>,
}

fn invalid_reference_data(
    reference: &FastaInput,
    record_name: &str,
    detail: impl std::fmt::Display,
) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidData,
        format!(
            "reference {} record {:?}: {detail}",
            reference.label, record_name
        ),
    )
}

fn extract_clean_sequence_minimizers(
    reference: &FastaInput,
    record_name: &str,
    sequence: &[u8],
    sequence_offset: usize,
    params: SketchParams,
) -> io::Result<Vec<ReferenceMinimizer>> {
    let SketchParams {
        kmer_size,
        window_size,
        minimizer_hash_seed,
        ..
    } = params;
    let Some(window_span) = kmer_size.checked_add(window_size.saturating_sub(1)) else {
        return Err(invalid_reference_data(
            reference,
            record_name,
            "minimizer window span exceeds usize",
        ));
    };
    if sequence.len() < window_span {
        return Ok(Vec::new());
    }

    let window_count: usize = sequence.len() - window_span + 1;
    if rayon::current_thread_index().is_none() || window_count <= EXTRACT_SEQUENCE_CHUNK_WINDOWS {
        return canonical_minimizers_with_positions(
            sequence,
            kmer_size,
            window_size,
            minimizer_hash_seed,
        )
        .into_iter()
        .map(|(hash, position)| {
            let global_position: usize = sequence_offset
                .checked_add(position as usize)
                .ok_or_else(|| {
                    invalid_reference_data(
                        reference,
                        record_name,
                        "minimizer position exceeds usize",
                    )
                })?;
            let position: u32 = u32::try_from(global_position).map_err(|err| {
                invalid_reference_data(
                    reference,
                    record_name,
                    format!(
                        "minimizer position {global_position} exceeds the u32 sketch limit: {err}"
                    ),
                )
            })?;
            Ok(ReferenceMinimizer { hash, position })
        })
        .collect();
    }

    // Partition by minimizer-window starts. Each noninitial chunk includes the
    // preceding window so adjacent equal minimizers are deduplicated exactly as
    // in a monolithic run; super-k-mer starts then assign every output once.
    let chunk_starts: Vec<usize> = (0..window_count)
        .step_by(EXTRACT_SEQUENCE_CHUNK_WINDOWS)
        .collect();
    let chunks: Vec<io::Result<Vec<ReferenceMinimizer>>> = chunk_starts
        .into_par_iter()
        .map(|core_start| {
            let core_end: usize = core_start
                .saturating_add(EXTRACT_SEQUENCE_CHUNK_WINDOWS)
                .min(window_count);
            let context_start: usize = core_start.saturating_sub(1);
            let context_end: usize = core_end.checked_add(window_span - 1).ok_or_else(|| {
                invalid_reference_data(reference, record_name, "chunk boundary exceeds usize")
            })?;
            let context: &[u8] = &sequence[context_start..context_end];
            canonical_minimizers_with_super_kmers(
                context,
                kmer_size,
                window_size,
                minimizer_hash_seed,
            )
            .into_iter()
            .filter_map(|(hash, position, super_kmer_start)| {
                let global_window_start: usize =
                    context_start + super_kmer_start as usize;
                (global_window_start >= core_start && global_window_start < core_end)
                    .then_some((hash, position))
            })
            .map(|(hash, position)| {
                let global_position: usize = sequence_offset
                    .checked_add(context_start)
                    .and_then(|offset| offset.checked_add(position as usize))
                    .ok_or_else(|| {
                        invalid_reference_data(
                            reference,
                            record_name,
                            "minimizer position exceeds usize",
                        )
                    })?;
                let position: u32 = u32::try_from(global_position).map_err(|err| {
                    invalid_reference_data(
                        reference,
                        record_name,
                        format!(
                            "minimizer position {global_position} exceeds the u32 sketch limit: {err}"
                        ),
                    )
                })?;
                Ok(ReferenceMinimizer { hash, position })
            })
            .collect()
        })
        .collect();
    let mut minimizers: Vec<ReferenceMinimizer> = Vec::new();
    for chunk in chunks {
        minimizers.extend(chunk?);
    }
    Ok(minimizers)
}

fn is_clean_base(base: u8) -> bool {
    matches!(base.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')
}

fn extract_segment_minimizers(
    reference: &FastaInput,
    record_name: &str,
    sequence: &[u8],
    params: SketchParams,
) -> io::Result<Vec<ReferenceMinimizer>> {
    let SketchParams {
        kmer_size,
        window_size,
        ..
    } = params;
    let Some(window_span) = kmer_size.checked_add(window_size.saturating_sub(1)) else {
        return Err(invalid_reference_data(
            reference,
            record_name,
            "minimizer window span exceeds usize",
        ));
    };
    if sequence.len() < window_span {
        return Ok(Vec::new());
    }

    let window_count: usize = sequence.len() - window_span + 1;
    let has_ambiguous_base: bool = sequence.iter().copied().any(|base| !is_clean_base(base));
    if !has_ambiguous_base {
        return extract_clean_sequence_minimizers(reference, record_name, sequence, 0, params);
    }
    if rayon::current_thread_index().is_none() || window_count <= EXTRACT_SEQUENCE_CHUNK_WINDOWS {
        return Ok(canonical_minimizers_with_positions(
            sequence,
            kmer_size,
            window_size,
            params.minimizer_hash_seed,
        )
        .into_iter()
        .map(|(hash, position)| ReferenceMinimizer { hash, position })
        .collect());
    }

    // An ambiguous base makes every minimizer window spanning it unusable. The library's
    // skip-ambiguous iterator also resets adjacent-value deduplication across that gap, so
    // independently extracting maximal clean ranges preserves the exact monolithic stream.
    let mut minimizers: Vec<ReferenceMinimizer> = Vec::new();
    let mut scan_position: usize = 0;
    while scan_position < sequence.len() {
        let mut ranges: Vec<Range<usize>> = Vec::with_capacity(EXTRACT_CLEAN_RANGE_BATCH);
        while scan_position < sequence.len() && ranges.len() < EXTRACT_CLEAN_RANGE_BATCH {
            while scan_position < sequence.len() && !is_clean_base(sequence[scan_position]) {
                scan_position += 1;
            }
            let clean_start: usize = scan_position;
            while scan_position < sequence.len() && is_clean_base(sequence[scan_position]) {
                scan_position += 1;
            }
            if scan_position.saturating_sub(clean_start) >= window_span {
                ranges.push(clean_start..scan_position);
            }
        }

        let range_results: Vec<io::Result<Vec<ReferenceMinimizer>>> = ranges
            .into_par_iter()
            .map(|range| {
                extract_clean_sequence_minimizers(
                    reference,
                    record_name,
                    &sequence[range.clone()],
                    range.start,
                    params,
                )
            })
            .collect();
        for range_result in range_results {
            minimizers.extend(range_result?);
        }
    }
    Ok(minimizers)
}

fn extract_record(
    reference: &FastaInput,
    record: PendingRecord,
    params: SketchParams,
) -> io::Result<ExtractedRecord> {
    let SketchParams {
        fragment_length,
        min_fragment_length,
        split_n_run,
        ..
    } = params;
    let mut mapped_length: u64 = 0;
    let mut segments: Vec<ExtractedReferenceSegment> = Vec::new();

    for segment_range in split_sequence_ranges(&record.sequence, split_n_run) {
        let segment_start: u32 = u32::try_from(segment_range.start).map_err(|err| {
            invalid_reference_data(
                reference,
                &record.name,
                format!(
                    "segment start {} exceeds the u32 sketch limit: {err}",
                    segment_range.start
                ),
            )
        })?;
        let segment_end: u32 = u32::try_from(segment_range.end).map_err(|err| {
            invalid_reference_data(
                reference,
                &record.name,
                format!(
                    "segment end {} exceeds the u32 sketch limit: {err}",
                    segment_range.end
                ),
            )
        })?;
        let segment_sequence: &[u8] = &record.sequence[segment_range];
        mapped_length = mapped_length
            .checked_add(mapped_length_from_fragment_ranges(
                segment_sequence.len(),
                fragment_length,
                min_fragment_length,
            ))
            .ok_or_else(|| {
                invalid_reference_data(reference, &record.name, "mapped length exceeds u64")
            })?;

        let mut minimizers: Vec<ReferenceMinimizer> =
            extract_segment_minimizers(reference, &record.name, segment_sequence, params)?;
        minimizers.sort_unstable_by_key(|minimizer| minimizer.position);

        segments.push(ExtractedReferenceSegment {
            record_name: record.name.clone(),
            segment_start,
            segment_end,
            minimizers,
        });
    }

    Ok(ExtractedRecord {
        mapped_length,
        segments,
    })
}

/// Read one FASTA input in bounded batches, extract each batch in the current Rayon pool, and
/// deliver segments to `sink` in their original record and split-segment order.
pub(crate) fn for_each_extracted_reference_segment<F>(
    reference: &FastaInput,
    params: SketchParams,
    mut sink: F,
) -> io::Result<ReferenceExtractionStats>
where
    F: FnMut(ExtractedReferenceSegment) -> io::Result<()>,
{
    let mut reader: fasta::io::Reader<Box<dyn io::BufRead>> = open_fasta_reader(&reference.open)
        .map_err(|err| {
            io::Error::new(
                err.kind(),
                format!("failed to open reference {}: {err}", reference.label),
            )
        })?;
    let mut records = reader.records();
    let mut record_number: usize = 0;
    let mut mapped_length: u64 = 0;
    let mut original_length: u64 = 0;

    loop {
        let mut batch: Vec<PendingRecord> = Vec::with_capacity(EXTRACT_BATCH_RECORDS);
        let mut batch_bases: usize = 0;

        while batch.len() < EXTRACT_BATCH_RECORDS && batch_bases < EXTRACT_BATCH_BASES {
            let Some(result) = records.next() else {
                break;
            };
            record_number = record_number.checked_add(1).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("reference {} has too many FASTA records", reference.label),
                )
            })?;
            let record: fasta::Record = result.map_err(|err| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "failed to read FASTA record {record_number} from reference {}: {err}",
                        reference.label
                    ),
                )
            })?;
            let sequence: Vec<u8> = record.sequence().as_ref().to_vec();
            original_length = original_length
                .checked_add(u64::try_from(sequence.len()).map_err(|err| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!(
                            "reference {} record length exceeds u64: {err}",
                            reference.label
                        ),
                    )
                })?)
                .ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("reference {} original length exceeds u64", reference.label),
                    )
                })?;
            batch_bases = batch_bases.checked_add(sequence.len()).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("reference {} batch size exceeds usize", reference.label),
                )
            })?;
            batch.push(PendingRecord {
                name: String::from_utf8_lossy(record.name()).into_owned(),
                sequence,
            });
        }

        if batch.is_empty() {
            break;
        }

        // Indexed parallel collection preserves the input order. When no Rayon worker is active
        // (the single-threaded direct-build path), avoid using the global pool.
        let extracted: Vec<io::Result<ExtractedRecord>> =
            if rayon::current_thread_index().is_some() && batch.len() > 1 {
                batch
                    .into_par_iter()
                    .map(|record| extract_record(reference, record, params))
                    .collect()
            } else {
                batch
                    .into_iter()
                    .map(|record| extract_record(reference, record, params))
                    .collect()
            };

        for extracted_record in extracted {
            let extracted_record: ExtractedRecord = extracted_record?;
            mapped_length = mapped_length
                .checked_add(extracted_record.mapped_length)
                .ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("reference {} mapped length exceeds u64", reference.label),
                    )
                })?;
            for segment in extracted_record.segments {
                sink(segment)?;
            }
        }
    }

    Ok(ReferenceExtractionStats {
        mapped_length,
        original_length,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ani::{
        DEFAULT_FRAGMENT_LENGTH, DEFAULT_KMER_SIZE, DEFAULT_MINIMIZER_HASH_SEED,
        DEFAULT_MIN_FRAGMENT_LENGTH, DEFAULT_WINDOW_SIZE,
    };
    use std::{env, fs, path::PathBuf, time::SystemTime};

    #[test]
    fn chunked_single_contig_extraction_matches_monolithic_output() -> io::Result<()> {
        let mut state: u64 = 0x9E37_79B9_7F4A_7C15;
        let sequence: Vec<u8> = (0..EXTRACT_SEQUENCE_CHUNK_WINDOWS + 4096)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state & 3) as usize]
            })
            .collect();
        let reference = FastaInput::from_path("chunked-reference.fa".to_owned());
        let params = SketchParams {
            kmer_size: DEFAULT_KMER_SIZE,
            window_size: DEFAULT_WINDOW_SIZE,
            minimizer_hash_seed: DEFAULT_MINIMIZER_HASH_SEED,
            fragment_length: DEFAULT_FRAGMENT_LENGTH,
            min_fragment_length: DEFAULT_MIN_FRAGMENT_LENGTH,
            split_n_run: 0,
        };
        let expected: Vec<ReferenceMinimizer> = canonical_minimizers_with_positions(
            &sequence,
            params.kmer_size,
            params.window_size,
            params.minimizer_hash_seed,
        )
        .into_iter()
        .map(|(hash, position)| ReferenceMinimizer { hash, position })
        .collect();
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .expect("test pool");
        let actual =
            pool.install(|| extract_segment_minimizers(&reference, "contig", &sequence, params))?;

        assert_eq!(actual, expected);

        let periodic: Vec<u8> = (0..EXTRACT_SEQUENCE_CHUNK_WINDOWS + 4096)
            .map(|position| b"ACGT"[position % 4])
            .collect();
        let periodic_expected: Vec<ReferenceMinimizer> = canonical_minimizers_with_positions(
            &periodic,
            params.kmer_size,
            params.window_size,
            params.minimizer_hash_seed,
        )
        .into_iter()
        .map(|(hash, position)| ReferenceMinimizer { hash, position })
        .collect();
        let periodic_actual = pool.install(|| {
            extract_segment_minimizers(&reference, "periodic-contig", &periodic, params)
        })?;
        assert_eq!(periodic_actual, periodic_expected);

        let mut ambiguous: Vec<u8> = periodic;
        // Exercise both a repeated minimizer on opposite sides of a skipped window and clean
        // ranges that themselves cross the chunking threshold.
        ambiguous[37] = b'N';
        ambiguous[99] = b'R';
        let ambiguous_expected: Vec<ReferenceMinimizer> = canonical_minimizers_with_positions(
            &ambiguous,
            params.kmer_size,
            params.window_size,
            params.minimizer_hash_seed,
        )
        .into_iter()
        .map(|(hash, position)| ReferenceMinimizer { hash, position })
        .collect();
        let ambiguous_actual = pool.install(|| {
            extract_segment_minimizers(&reference, "ambiguous-contig", &ambiguous, params)
        })?;
        assert_eq!(ambiguous_actual, ambiguous_expected);
        Ok(())
    }

    #[test]
    fn extraction_preserves_order_across_batches_and_split_segments() -> io::Result<()> {
        let reference_path: PathBuf = env::temp_dir().join(format!(
            "fasterani_extract_order_{}_{}.fa",
            std::process::id(),
            SystemTime::now()
                .duration_since(SystemTime::UNIX_EPOCH)
                .expect("clock after epoch")
                .as_nanos()
        ));
        let mut fasta = String::new();
        for record_index in 0..(EXTRACT_BATCH_RECORDS + 3) {
            fasta.push_str(&format!(">record-{record_index}\nACGTACGTNNNNACGTACGT\n"));
        }
        fs::write(&reference_path, fasta)?;
        let reference = FastaInput::from_path(reference_path.to_string_lossy().into_owned());
        let params = SketchParams {
            kmer_size: DEFAULT_KMER_SIZE,
            window_size: DEFAULT_WINDOW_SIZE,
            minimizer_hash_seed: DEFAULT_MINIMIZER_HASH_SEED,
            fragment_length: DEFAULT_FRAGMENT_LENGTH,
            min_fragment_length: DEFAULT_MIN_FRAGMENT_LENGTH,
            split_n_run: 4,
        };
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .expect("test pool");
        let mut observed = Vec::new();
        let extraction = pool.install(|| {
            for_each_extracted_reference_segment(&reference, params, |segment| {
                observed.push((
                    segment.record_name,
                    segment.segment_start,
                    segment.segment_end,
                ));
                Ok(())
            })
        })?;

        assert_eq!(extraction.mapped_length, 0);
        assert_eq!(
            extraction.original_length,
            ((EXTRACT_BATCH_RECORDS + 3) * 20) as u64
        );
        assert_eq!(observed.len(), (EXTRACT_BATCH_RECORDS + 3) * 2);
        for record_index in 0..(EXTRACT_BATCH_RECORDS + 3) {
            assert_eq!(
                observed[record_index * 2],
                (format!("record-{record_index}"), 0, 8)
            );
            assert_eq!(
                observed[record_index * 2 + 1],
                (format!("record-{record_index}"), 12, 20)
            );
        }

        fs::remove_file(reference_path)?;
        Ok(())
    }
}
