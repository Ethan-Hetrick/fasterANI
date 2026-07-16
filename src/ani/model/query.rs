//! Query-side data types: fragments, mapping results, and ANI summaries.

use std::collections::HashSet;
use std::io::{self, Read, Write};
use std::mem::size_of;

use crate::ani::{MinimizerKey, SeedHit, SlidingSketchCounter};

/// One fixed-length query fragment with full scoring minimizers and optional reduced seed hashes.
#[derive(Clone)]
pub(crate) struct QueryFragment {
    pub(crate) id: usize,
    pub(crate) contig_id: usize,
    pub(crate) start: u32,
    pub(crate) end: u32,
    pub(crate) length: u32,
    pub(crate) minimizers: Vec<MinimizerKey>,
    pub(crate) seed_minimizers: Vec<MinimizerKey>,
}

/// Query genome split into FastANI-style fixed-width fragments.
pub(crate) struct QueryFile {
    pub(crate) fragments: Vec<QueryFragment>,
    pub(crate) contig_names: Vec<String>,
    pub(crate) mapped_length: u64,
}

const QUERY_SPOOL_MAGIC: [u8; 8] = *b"FANIQRY\0";
const QUERY_SPOOL_VERSION: u32 = 1;
// Query spools are internal temporary files, but keep corrupt inputs from requesting
// effectively unbounded allocations before an EOF or format error can be reported.
const MAX_QUERY_SPOOL_ALLOCATION_BYTES: u64 = 16 * 1024 * 1024 * 1024;

impl QueryFile {
    /// Write a deterministic, architecture-independent representation of this prepared query.
    pub(crate) fn write_to<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        validate_spool_allocation_budget(self)?;

        write_all(writer, &QUERY_SPOOL_MAGIC, "magic")?;
        write_u32(writer, QUERY_SPOOL_VERSION, "format version")?;
        write_u64(writer, self.mapped_length, "mapped length")?;
        write_len(writer, self.contig_names.len(), "contig count")?;
        write_len(writer, self.fragments.len(), "fragment count")?;

        for (contig_index, name) in self.contig_names.iter().enumerate() {
            write_len(
                writer,
                name.len(),
                &format!("contig {contig_index} name length"),
            )?;
            write_all(
                writer,
                name.as_bytes(),
                &format!("contig {contig_index} name"),
            )?;
        }

        for (fragment_index, fragment) in self.fragments.iter().enumerate() {
            write_usize_as_u64(
                writer,
                fragment.id,
                &format!("fragment {fragment_index} id"),
            )?;
            write_usize_as_u64(
                writer,
                fragment.contig_id,
                &format!("fragment {fragment_index} contig id"),
            )?;
            write_u32(
                writer,
                fragment.start,
                &format!("fragment {fragment_index} start"),
            )?;
            write_u32(
                writer,
                fragment.end,
                &format!("fragment {fragment_index} end"),
            )?;
            write_u32(
                writer,
                fragment.length,
                &format!("fragment {fragment_index} length"),
            )?;
            write_minimizers(
                writer,
                &fragment.minimizers,
                &format!("fragment {fragment_index} minimizers"),
            )?;
            write_minimizers(
                writer,
                &fragment.seed_minimizers,
                &format!("fragment {fragment_index} seed minimizers"),
            )?;
        }

        Ok(())
    }

    /// Read one complete prepared query spool, rejecting unsupported or malformed data.
    pub(crate) fn read_from<R: Read>(reader: &mut R) -> io::Result<Self> {
        let mut magic = [0_u8; QUERY_SPOOL_MAGIC.len()];
        read_exact(reader, &mut magic, "magic")?;
        if magic != QUERY_SPOOL_MAGIC {
            return Err(invalid_data("query spool has an invalid magic header"));
        }

        let version = read_u32(reader, "format version")?;
        if version != QUERY_SPOOL_VERSION {
            return Err(invalid_data(format!(
                "unsupported query spool version {version}; expected {QUERY_SPOOL_VERSION}"
            )));
        }

        let mapped_length = read_u64(reader, "mapped length")?;
        let contig_count = read_len(reader, "contig count")?;
        let fragment_count = read_len(reader, "fragment count")?;
        let mut budget = AllocationBudget::new();
        budget.charge::<String>(contig_count.encoded, "contig name table")?;
        budget.charge::<QueryFragment>(fragment_count.encoded, "fragment table")?;

        let mut contig_names = Vec::new();
        reserve_exact(&mut contig_names, contig_count.decoded, "contig name table")?;
        for contig_index in 0..contig_count.decoded {
            let context = format!("contig {contig_index} name");
            let name_length = read_len(reader, &format!("{context} length"))?;
            budget.charge::<u8>(name_length.encoded, &context)?;
            let mut bytes = Vec::new();
            reserve_exact(&mut bytes, name_length.decoded, &context)?;
            bytes.resize(name_length.decoded, 0);
            read_exact(reader, &mut bytes, &context)?;
            let name = String::from_utf8(bytes).map_err(|error| {
                invalid_data(format!("query spool {context} is not valid UTF-8: {error}"))
            })?;
            contig_names.push(name);
        }

        let mut fragments = Vec::new();
        reserve_exact(&mut fragments, fragment_count.decoded, "fragment table")?;
        for fragment_index in 0..fragment_count.decoded {
            let prefix = format!("fragment {fragment_index}");
            let id = read_usize(reader, &format!("{prefix} id"))?;
            let contig_id = read_usize(reader, &format!("{prefix} contig id"))?;
            if contig_id >= contig_names.len() {
                return Err(invalid_data(format!(
                    "query spool {prefix} contig id {contig_id} is out of range for {} contigs",
                    contig_names.len()
                )));
            }
            let start = read_u32(reader, &format!("{prefix} start"))?;
            let end = read_u32(reader, &format!("{prefix} end"))?;
            let length = read_u32(reader, &format!("{prefix} length"))?;
            let minimizers = read_minimizers(reader, &mut budget, &format!("{prefix} minimizers"))?;
            let seed_minimizers =
                read_minimizers(reader, &mut budget, &format!("{prefix} seed minimizers"))?;
            fragments.push(QueryFragment {
                id,
                contig_id,
                start,
                end,
                length,
                minimizers,
                seed_minimizers,
            });
        }

        reject_trailing_bytes(reader)?;
        Ok(Self {
            fragments,
            contig_names,
            mapped_length,
        })
    }
}

struct EncodedLength {
    encoded: u64,
    decoded: usize,
}

struct AllocationBudget {
    remaining: u64,
}

impl AllocationBudget {
    fn new() -> Self {
        Self {
            remaining: MAX_QUERY_SPOOL_ALLOCATION_BYTES,
        }
    }

    fn charge<T>(&mut self, count: u64, context: &str) -> io::Result<()> {
        let element_size = u64::try_from(size_of::<T>()).map_err(|_| {
            invalid_data(format!(
                "query spool {context} element size does not fit in u64"
            ))
        })?;
        let bytes = count.checked_mul(element_size).ok_or_else(|| {
            invalid_data(format!("query spool {context} allocation size overflow"))
        })?;
        self.remaining = self.remaining.checked_sub(bytes).ok_or_else(|| {
            invalid_data(format!(
                "query spool allocations exceed the {}-byte safety limit while reading {context}",
                MAX_QUERY_SPOOL_ALLOCATION_BYTES
            ))
        })?;
        Ok(())
    }
}

fn validate_spool_allocation_budget(query: &QueryFile) -> io::Result<()> {
    let mut budget = AllocationBudget::new();
    budget.charge::<String>(
        usize_to_u64(query.contig_names.len(), "contig count")?,
        "contig name table",
    )?;
    budget.charge::<QueryFragment>(
        usize_to_u64(query.fragments.len(), "fragment count")?,
        "fragment table",
    )?;
    for (contig_index, name) in query.contig_names.iter().enumerate() {
        budget.charge::<u8>(
            usize_to_u64(name.len(), &format!("contig {contig_index} name length"))?,
            &format!("contig {contig_index} name"),
        )?;
    }
    for (fragment_index, fragment) in query.fragments.iter().enumerate() {
        budget.charge::<MinimizerKey>(
            usize_to_u64(
                fragment.minimizers.len(),
                &format!("fragment {fragment_index} minimizer count"),
            )?,
            &format!("fragment {fragment_index} minimizers"),
        )?;
        budget.charge::<MinimizerKey>(
            usize_to_u64(
                fragment.seed_minimizers.len(),
                &format!("fragment {fragment_index} seed minimizer count"),
            )?,
            &format!("fragment {fragment_index} seed minimizers"),
        )?;
    }
    Ok(())
}

fn write_minimizers<W: Write>(
    writer: &mut W,
    minimizers: &[MinimizerKey],
    context: &str,
) -> io::Result<()> {
    write_len(writer, minimizers.len(), &format!("{context} count"))?;
    for minimizer in minimizers.iter().copied() {
        write_u32(writer, minimizer, context)?;
    }
    Ok(())
}

fn read_minimizers<R: Read>(
    reader: &mut R,
    budget: &mut AllocationBudget,
    context: &str,
) -> io::Result<Vec<MinimizerKey>> {
    let count = read_len(reader, &format!("{context} count"))?;
    budget.charge::<MinimizerKey>(count.encoded, context)?;
    let mut minimizers = Vec::new();
    reserve_exact(&mut minimizers, count.decoded, context)?;
    for _ in 0..count.decoded {
        minimizers.push(read_u32(reader, context)?);
    }
    Ok(minimizers)
}

fn reserve_exact<T>(values: &mut Vec<T>, count: usize, context: &str) -> io::Result<()> {
    values.try_reserve_exact(count).map_err(|error| {
        invalid_data(format!(
            "cannot allocate query spool {context} with {count} items: {error}"
        ))
    })
}

fn write_len<W: Write>(writer: &mut W, value: usize, context: &str) -> io::Result<()> {
    write_usize_as_u64(writer, value, context)
}

fn write_usize_as_u64<W: Write>(writer: &mut W, value: usize, context: &str) -> io::Result<()> {
    write_u64(writer, usize_to_u64(value, context)?, context)
}

fn usize_to_u64(value: usize, context: &str) -> io::Result<u64> {
    u64::try_from(value).map_err(|_| {
        invalid_data(format!(
            "query spool {context} value {value} does not fit in u64"
        ))
    })
}

fn read_len<R: Read>(reader: &mut R, context: &str) -> io::Result<EncodedLength> {
    let encoded = read_u64(reader, context)?;
    let decoded = usize::try_from(encoded).map_err(|_| {
        invalid_data(format!(
            "query spool {context} value {encoded} does not fit in usize"
        ))
    })?;
    Ok(EncodedLength { encoded, decoded })
}

fn read_usize<R: Read>(reader: &mut R, context: &str) -> io::Result<usize> {
    Ok(read_len(reader, context)?.decoded)
}

fn write_u32<W: Write>(writer: &mut W, value: u32, context: &str) -> io::Result<()> {
    write_all(writer, &value.to_le_bytes(), context)
}

fn write_u64<W: Write>(writer: &mut W, value: u64, context: &str) -> io::Result<()> {
    write_all(writer, &value.to_le_bytes(), context)
}

fn read_u32<R: Read>(reader: &mut R, context: &str) -> io::Result<u32> {
    let mut bytes = [0_u8; size_of::<u32>()];
    read_exact(reader, &mut bytes, context)?;
    Ok(u32::from_le_bytes(bytes))
}

fn read_u64<R: Read>(reader: &mut R, context: &str) -> io::Result<u64> {
    let mut bytes = [0_u8; size_of::<u64>()];
    read_exact(reader, &mut bytes, context)?;
    Ok(u64::from_le_bytes(bytes))
}

fn write_all<W: Write>(writer: &mut W, bytes: &[u8], context: &str) -> io::Result<()> {
    writer.write_all(bytes).map_err(|error| {
        io::Error::new(
            error.kind(),
            format!("failed to write query spool {context}: {error}"),
        )
    })
}

fn read_exact<R: Read>(reader: &mut R, bytes: &mut [u8], context: &str) -> io::Result<()> {
    reader.read_exact(bytes).map_err(|error| {
        let kind = if error.kind() == io::ErrorKind::UnexpectedEof {
            io::ErrorKind::UnexpectedEof
        } else {
            error.kind()
        };
        io::Error::new(
            kind,
            format!("failed to read query spool {context}: {error}"),
        )
    })
}

fn reject_trailing_bytes<R: Read>(reader: &mut R) -> io::Result<()> {
    let mut byte = [0_u8; 1];
    loop {
        match reader.read(&mut byte) {
            Ok(0) => return Ok(()),
            Ok(_) => return Err(invalid_data("query spool contains trailing bytes")),
            Err(error) if error.kind() == io::ErrorKind::Interrupted => continue,
            Err(error) => {
                return Err(io::Error::new(
                    error.kind(),
                    format!("failed to check the end of the query spool: {error}"),
                ));
            }
        }
    }
}

fn invalid_data(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.into())
}

/// Candidate reference position interval produced by clustered seed hits.
#[derive(Clone, Copy)]
pub(crate) struct ReferenceCandidateRegion {
    pub(crate) reference_contig_id: usize,
    pub(crate) start_position: u32,
    pub(crate) end_position: u32,
}

/// One retained mapping between a query fragment and a reference location.
#[derive(Clone)]
pub(crate) struct MappingResult {
    pub(crate) reference_file_id: usize,
    pub(crate) reference_contig_id: usize,
    pub(crate) query_fragment_id: usize,
    pub(crate) query_fragment_length: u32,
    pub(crate) reference_start: u32,
    pub(crate) identity: f64,
    pub(crate) query_minimizer_count: usize,
    pub(crate) reference_minimizer_count: usize,
    pub(crate) shared_minimizers: usize,
    pub(crate) union_minimizers: usize,
    pub(crate) jaccard: f64,
}

/// Stable identity for marking mappings that survive reciprocal-best filtering.
#[derive(Hash, Eq, PartialEq)]
pub(crate) struct MappingResultKey {
    pub(crate) reference_file_id: usize,
    pub(crate) reference_contig_id: usize,
    pub(crate) query_fragment_id: usize,
    pub(crate) query_fragment_length: u32,
    pub(crate) reference_start: u32,
    pub(crate) identity_bits: u64,
    pub(crate) query_minimizer_count: usize,
    pub(crate) reference_minimizer_count: usize,
    pub(crate) shared_minimizers: usize,
    pub(crate) union_minimizers: usize,
}

impl MappingResultKey {
    pub(crate) fn from_mapping(mapping: &MappingResult) -> Self {
        Self {
            reference_file_id: mapping.reference_file_id,
            reference_contig_id: mapping.reference_contig_id,
            query_fragment_id: mapping.query_fragment_id,
            query_fragment_length: mapping.query_fragment_length,
            reference_start: mapping.reference_start,
            identity_bits: mapping.identity.to_bits(),
            query_minimizer_count: mapping.query_minimizer_count,
            reference_minimizer_count: mapping.reference_minimizer_count,
            shared_minimizers: mapping.shared_minimizers,
            union_minimizers: mapping.union_minimizers,
        }
    }
}

/// Accumulator for final ANI output for one reference file.
#[derive(Clone, Default)]
pub(crate) struct AniSummary {
    pub(crate) shared_fragments: usize,
    pub(crate) shared_bases: u64,
    pub(crate) weighted_identity_sum: f64,
    pub(crate) distribution_stats: AniDistributionStats,
}

/// Final ANI output for one query contig against one reference file.
#[derive(Clone, Default)]
pub(crate) struct ContigAniSummary {
    pub(crate) eligible_fragments: usize,
    pub(crate) summary: AniSummary,
}

/// Distribution statistics for retained fragment ANI values.
#[derive(Clone, Copy)]
pub struct AniDistributionStats {
    pub median: f64,
    pub stddev: f64,
    pub mad: f64,
    pub ci_95_lower: f64,
    pub ci_95_upper: f64,
    pub f99: f64,
    // P99/P80 are kept for future evaluation, but are intentionally not part
    // of the public TSV output until their interpretation is settled.
    #[allow(dead_code)]
    pub p99: f64,
    pub f80: f64,
    #[allow(dead_code)]
    pub p80: f64,
}

impl Default for AniDistributionStats {
    fn default() -> Self {
        Self {
            median: f64::NAN,
            stddev: f64::NAN,
            mad: f64::NAN,
            ci_95_lower: f64::NAN,
            ci_95_upper: f64::NAN,
            f99: f64::NAN,
            p99: f64::NAN,
            f80: f64::NAN,
            p80: f64::NAN,
        }
    }
}

/// Final per-reference summaries plus the exact mappings that contributed to them.
pub(crate) struct AniComputation {
    pub(crate) summaries: Vec<AniSummary>,
    pub(crate) contig_summaries: Vec<Vec<ContigAniSummary>>,
    pub(crate) reciprocal_best_keys: HashSet<MappingResultKey>,
}

/// Per-thread reusable buffers for query-to-reference mapping.
#[derive(Default)]
pub(crate) struct MappingScratch {
    pub(crate) seed_hits: Vec<SeedHit>,
    pub(crate) candidate_regions: Vec<ReferenceCandidateRegion>,
    pub(crate) fragment_mappings: Vec<MappingResult>,
    pub(crate) counter: SlidingSketchCounter,
    pub(crate) slot_sorted_minimizers: Vec<(u64, MinimizerKey)>,
    pub(crate) hit_ranges: Vec<(usize, usize)>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn sample_query() -> QueryFile {
        QueryFile {
            fragments: vec![
                QueryFragment {
                    id: 7,
                    contig_id: 1,
                    start: 11,
                    end: 29,
                    length: 18,
                    minimizers: vec![u32::MAX, 0, 42],
                    seed_minimizers: vec![9, 3],
                },
                QueryFragment {
                    id: 2,
                    contig_id: 0,
                    start: 0,
                    end: 5,
                    length: 5,
                    minimizers: Vec::new(),
                    seed_minimizers: vec![17],
                },
            ],
            contig_names: vec!["first contig".to_owned(), "β-contig".to_owned()],
            mapped_length: u64::from(u32::MAX) + 19,
        }
    }

    fn encode(query: &QueryFile) -> Vec<u8> {
        let mut bytes = Vec::new();
        query.write_to(&mut bytes).unwrap();
        bytes
    }

    fn assert_query_eq(actual: &QueryFile, expected: &QueryFile) {
        assert_eq!(actual.mapped_length, expected.mapped_length);
        assert_eq!(actual.contig_names, expected.contig_names);
        assert_eq!(actual.fragments.len(), expected.fragments.len());
        for (actual, expected) in actual.fragments.iter().zip(&expected.fragments) {
            assert_eq!(actual.id, expected.id);
            assert_eq!(actual.contig_id, expected.contig_id);
            assert_eq!(actual.start, expected.start);
            assert_eq!(actual.end, expected.end);
            assert_eq!(actual.length, expected.length);
            assert_eq!(actual.minimizers, expected.minimizers);
            assert_eq!(actual.seed_minimizers, expected.seed_minimizers);
        }
    }

    #[test]
    fn query_spool_round_trip_preserves_every_field_and_order() {
        let expected = sample_query();
        let bytes = encode(&expected);
        let actual = QueryFile::read_from(&mut Cursor::new(bytes)).unwrap();
        assert_query_eq(&actual, &expected);
    }

    #[test]
    fn query_spool_encoding_is_deterministic_and_little_endian() {
        let query = sample_query();
        let first = encode(&query);
        let second = encode(&query);

        assert_eq!(first, second);
        assert_eq!(&first[..8], &QUERY_SPOOL_MAGIC);
        assert_eq!(&first[8..12], &QUERY_SPOOL_VERSION.to_le_bytes());
        assert_eq!(&first[12..20], &query.mapped_length.to_le_bytes());
        assert_eq!(&first[20..28], &2_u64.to_le_bytes());
        assert_eq!(&first[28..36], &2_u64.to_le_bytes());
    }

    #[test]
    fn query_spool_reports_truncated_field_context() {
        let mut bytes = encode(&sample_query());
        bytes.pop();

        let error = QueryFile::read_from(&mut Cursor::new(bytes))
            .err()
            .expect("truncated spool should fail");
        assert_eq!(error.kind(), io::ErrorKind::UnexpectedEof);
        assert!(error.to_string().contains("fragment 1 seed minimizers"));
    }

    #[test]
    fn query_spool_rejects_bad_magic_and_unsupported_version() {
        let mut bad_magic = encode(&sample_query());
        bad_magic[0] ^= 0xff;
        let error = QueryFile::read_from(&mut Cursor::new(bad_magic))
            .err()
            .expect("bad magic should fail");
        assert_eq!(error.kind(), io::ErrorKind::InvalidData);
        assert!(error.to_string().contains("magic"));

        let mut bad_version = encode(&sample_query());
        bad_version[8..12].copy_from_slice(&(QUERY_SPOOL_VERSION + 1).to_le_bytes());
        let error = QueryFile::read_from(&mut Cursor::new(bad_version))
            .err()
            .expect("unsupported version should fail");
        assert_eq!(error.kind(), io::ErrorKind::InvalidData);
        assert!(error
            .to_string()
            .contains("unsupported query spool version"));
    }

    #[test]
    fn query_spool_rejects_allocation_over_budget_before_reserving() {
        let string_size = u64::try_from(size_of::<String>()).unwrap();
        let excessive_contig_count = MAX_QUERY_SPOOL_ALLOCATION_BYTES / string_size + 1;
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&QUERY_SPOOL_MAGIC);
        bytes.extend_from_slice(&QUERY_SPOOL_VERSION.to_le_bytes());
        bytes.extend_from_slice(&0_u64.to_le_bytes());
        bytes.extend_from_slice(&excessive_contig_count.to_le_bytes());
        bytes.extend_from_slice(&0_u64.to_le_bytes());

        let error = QueryFile::read_from(&mut Cursor::new(bytes))
            .err()
            .expect("over-budget allocation should fail");
        assert_eq!(error.kind(), io::ErrorKind::InvalidData);
        assert!(error.to_string().contains("safety limit"));
    }

    #[test]
    fn query_spool_rejects_out_of_range_fragment_contig_id() {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&QUERY_SPOOL_MAGIC);
        bytes.extend_from_slice(&QUERY_SPOOL_VERSION.to_le_bytes());
        bytes.extend_from_slice(&0_u64.to_le_bytes());
        bytes.extend_from_slice(&0_u64.to_le_bytes());
        bytes.extend_from_slice(&1_u64.to_le_bytes());
        bytes.extend_from_slice(&0_u64.to_le_bytes());
        bytes.extend_from_slice(&0_u64.to_le_bytes());

        let error = QueryFile::read_from(&mut Cursor::new(bytes))
            .err()
            .expect("out-of-range contig id should fail");
        assert_eq!(error.kind(), io::ErrorKind::InvalidData);
        assert!(error.to_string().contains("contig id 0 is out of range"));
    }

    #[test]
    fn query_spool_rejects_invalid_utf8_and_trailing_bytes() {
        let mut invalid_utf8 = encode(&QueryFile {
            fragments: Vec::new(),
            contig_names: vec!["x".to_owned()],
            mapped_length: 0,
        });
        // Header is 36 bytes, followed by the first name length and then its bytes.
        invalid_utf8[44] = 0xff;
        let error = QueryFile::read_from(&mut Cursor::new(invalid_utf8))
            .err()
            .expect("invalid UTF-8 should fail");
        assert_eq!(error.kind(), io::ErrorKind::InvalidData);
        assert!(error.to_string().contains("UTF-8"));

        let mut trailing = encode(&sample_query());
        trailing.push(0);
        let error = QueryFile::read_from(&mut Cursor::new(trailing))
            .err()
            .expect("trailing bytes should fail");
        assert_eq!(error.kind(), io::ErrorKind::InvalidData);
        assert!(error.to_string().contains("trailing bytes"));
    }
}
