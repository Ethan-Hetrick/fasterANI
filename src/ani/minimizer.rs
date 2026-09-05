//! Canonical minimizer extraction, sliding-window sketches, and query fragmentation.

use std::io;

use seq_hash::NtHasher;
use simd_minimizers::canonical_minimizers;
use simd_minimizers::packed_seq::{PackedNSeqVec, Seq};

use crate::ani::{constants::MinimizerKey, model::reference::ReferenceMinimizer};

/// Sliding minimizer set used while scoring candidate reference windows.
#[derive(Default)]
pub(crate) struct SlidingSketchCounter {
    pub(crate) hashes: Vec<MinimizerKey>,
    pub(crate) query_present: Vec<bool>,
    pub(crate) reference_counts: Vec<usize>,
    pub(crate) active: Fenwick,
    pub(crate) shared: Fenwick,
    pub(crate) sketch_size: usize,
}

/// Minimizers observed for one sequence plus the raw minimizer-window coverage signal.
pub(crate) struct MinimizerObservation {
    pub(crate) minimizers_with_positions: Vec<(MinimizerKey, u32)>,
}

/// Query-fragment sketch with both scoring minimizers and the missing-window quality signal.
pub(crate) struct QueryFragmentSketch {
    pub(crate) minimizers: Vec<MinimizerKey>,
    pub(crate) seed_minimizers: Vec<MinimizerKey>,
}

/// Fenwick tree for rank/select operations over active minimizer coordinates.
#[derive(Default)]
pub(crate) struct Fenwick {
    pub(crate) tree: Vec<i32>,
    pub(crate) total: i32,
}

impl Fenwick {
    pub(crate) fn reset(&mut self, len: usize) {
        let tree_len: usize = len.saturating_add(1);
        self.tree.clear();
        self.tree.reserve(tree_len);
        self.tree.resize(tree_len, 0);
        self.total = 0;
    }

    pub(crate) fn add(&mut self, index: usize, delta: i32) {
        let mut i = index + 1;
        self.total += delta;

        while i < self.tree.len() {
            self.tree[i] += delta;
            i += i & i.wrapping_neg();
        }
    }

    pub(crate) fn prefix_sum_inclusive(&self, index: usize) -> i32 {
        let mut i = index + 1;
        let mut sum = 0;

        while i > 0 {
            sum += self.tree[i];
            i &= i - 1;
        }

        sum
    }

    pub(crate) fn select_by_rank(&self, rank: i32) -> Option<usize> {
        if rank <= 0 || rank > self.total {
            return None;
        }

        let mut idx = 0usize;
        let mut bit = 1usize;

        while bit < self.tree.len() {
            bit <<= 1;
        }

        let mut remaining = rank;

        while bit > 0 {
            let next = idx + bit;

            if next < self.tree.len() && self.tree[next] < remaining {
                idx = next;
                remaining -= self.tree[next];
            }

            bit >>= 1;
        }

        Some(idx)
    }
}

impl SlidingSketchCounter {
    pub(crate) fn clear(&mut self) {
        for count in &mut self.reference_counts {
            *count = 0;
        }

        self.active.reset(self.hashes.len());
        self.shared.reset(self.hashes.len());

        for (index, &present) in self.query_present.iter().enumerate() {
            if present {
                self.active.add(index, 1);
            }
        }
    }

    pub(crate) fn prepare(
        &mut self,
        query_minimizers: &[MinimizerKey],
        reference_minimizers: &[ReferenceMinimizer],
    ) {
        let coordinate_capacity: usize = query_minimizers
            .len()
            .saturating_add(reference_minimizers.len());

        self.hashes.clear();
        self.hashes.reserve(coordinate_capacity);
        self.hashes.extend_from_slice(query_minimizers);
        self.hashes
            .extend(reference_minimizers.iter().map(|minimizer| minimizer.hash));
        self.hashes.sort_unstable();
        self.hashes.dedup();

        let coordinate_count: usize = self.hashes.len();
        self.query_present.clear();
        self.query_present.reserve(coordinate_count);
        self.query_present.resize(coordinate_count, false);

        for minimizer in query_minimizers {
            let index = self
                .hashes
                .binary_search(minimizer)
                .expect("query minimizer missing from coordinate table");
            self.query_present[index] = true;
        }

        self.reference_counts.clear();
        self.reference_counts.reserve(coordinate_count);
        self.reference_counts.resize(coordinate_count, 0);
        self.sketch_size = query_minimizers.len();
        self.clear();
    }

    pub(crate) fn insert(&mut self, hash: MinimizerKey) {
        // Keep this as binary search unless default fragment geometry is re-profiled:
        // the reverted SIMD scan assumed much smaller coordinate tables than the
        // 150-300 entry common case and added enough branch-mispredict overhead
        // in this loop to outweigh the intended lookup savings.
        let index = self
            .hashes
            .binary_search(&hash)
            .expect("reference minimizer missing from coordinate table");
        let count = &mut self.reference_counts[index];

        if *count == 0 {
            if self.query_present[index] {
                self.shared.add(index, 1);
            } else {
                self.active.add(index, 1);
            }
        }

        *count += 1;
    }

    pub(crate) fn remove(&mut self, hash: MinimizerKey) {
        let index = self
            .hashes
            .binary_search(&hash)
            .expect("reference minimizer missing from coordinate table");
        let count = &mut self.reference_counts[index];

        if *count == 0 {
            return;
        }

        *count -= 1;

        if *count == 0 {
            if self.query_present[index] {
                self.shared.add(index, -1);
            } else {
                self.active.add(index, -1);
            }
        }
    }

    pub(crate) fn shared_count(&self) -> usize {
        let rank = self.sketch_size as i32;
        let Some(pivot) = self.active.select_by_rank(rank) else {
            return 0;
        };

        self.shared.prefix_sum_inclusive(pivot) as usize
    }

    pub(crate) fn reference_minimizer_count(&self) -> usize {
        self.reference_counts
            .iter()
            .filter(|&&count| count > 0)
            .count()
    }
}

pub(crate) fn expected_minimizer_window_count(
    sequence_len: usize,
    kmer_size: usize,
    window_size: usize,
) -> usize {
    let Some(minimum_sequence_len) = kmer_size.checked_add(window_size.saturating_sub(1)) else {
        return 0;
    };

    if sequence_len < minimum_sequence_len {
        0
    } else {
        sequence_len - minimum_sequence_len + 1
    }
}

#[cfg(test)]
pub(crate) fn usable_minimizer_window_count(
    sequence: &[u8],
    kmer_size: usize,
    window_size: usize,
) -> usize {
    let expected_window_count: usize =
        expected_minimizer_window_count(sequence.len(), kmer_size, window_size);
    if expected_window_count == 0 {
        return 0;
    }

    let minimizer_window_span: usize = kmer_size + window_size - 1;
    let mut ambiguous_prefix: Vec<usize> = Vec::with_capacity(sequence.len() + 1);
    ambiguous_prefix.push(0);

    for base in sequence {
        let previous_count: usize = *ambiguous_prefix
            .last()
            .expect("ambiguous prefix always has a zero entry");
        let next_count: usize = previous_count + usize::from(matches!(base, b'N' | b'n'));
        ambiguous_prefix.push(next_count);
    }

    (0..expected_window_count)
        .filter(|&start| ambiguous_prefix[start + minimizer_window_span] == ambiguous_prefix[start])
        .count()
}

/// Return canonical minimizers plus how many minimizer windows were usable before deduplication.
pub(crate) fn canonical_minimizer_observation(
    sequence: &[u8],
    kmer_size: usize,
    window_size: usize,
    minimizer_hash_seed: u32,
) -> MinimizerObservation {
    if sequence.len() < kmer_size + window_size.saturating_sub(1) || is_all_n_sequence(sequence) {
        return MinimizerObservation {
            minimizers_with_positions: Vec::new(),
        };
    }

    let packed_sequence: PackedNSeqVec = PackedNSeqVec::from_ascii(sequence);
    let packed_sequence_slice = packed_sequence.as_slice();
    let sequence_slice = packed_sequence_slice.seq;
    let hasher: NtHasher<true> = NtHasher::<true>::new_with_seed(kmer_size, minimizer_hash_seed);
    let mut minimizer_positions: Vec<u32> = Vec::new();
    let minimizer_builder = canonical_minimizers(kmer_size, window_size).hasher(&hasher);
    let _ = minimizer_builder
        .run_skip_ambiguous_windows(packed_sequence_slice, &mut minimizer_positions);
    debug_assert!(kmer_size <= 16, "u32 2-bit minimizer keys require k <= 16");

    let minimizers_with_positions: Vec<(MinimizerKey, u32)> = minimizer_positions
        .into_iter()
        .filter_map(|position| {
            let pos = position as usize;
            let forward_kmer = sequence_slice.read_kmer(kmer_size, pos);
            let reverse_complement_kmer = sequence_slice.read_revcomp_kmer(kmer_size, pos);
            if forward_kmer == reverse_complement_kmer {
                return None;
            }

            let canonical_kmer: u64 = forward_kmer.min(reverse_complement_kmer);
            let key: MinimizerKey = MinimizerKey::try_from(canonical_kmer)
                .expect("u32 2-bit minimizer keys require k <= 16");

            Some((key, position))
        })
        .collect();

    MinimizerObservation {
        minimizers_with_positions,
    }
}

/// Return canonical minimizer hashes and their positions for one nucleotide sequence.
pub(crate) fn canonical_minimizers_with_positions(
    sequence: &[u8],
    kmer_size: usize,
    window_size: usize,
    minimizer_hash_seed: u32,
) -> Vec<(MinimizerKey, u32)> {
    canonical_minimizer_observation(sequence, kmer_size, window_size, minimizer_hash_seed)
        .minimizers_with_positions
}

/// Return canonical minimizers together with the first window start at which
/// each deduplicated minimizer was selected. This clean-sequence variant is
/// used to split a long contig into independently computable window ranges.
pub(crate) fn canonical_minimizers_with_super_kmers(
    sequence: &[u8],
    kmer_size: usize,
    window_size: usize,
    minimizer_hash_seed: u32,
) -> Vec<(MinimizerKey, u32, u32)> {
    if sequence.len() < kmer_size + window_size.saturating_sub(1) {
        return Vec::new();
    }
    debug_assert!(
        !sequence.iter().any(|base| matches!(base, b'N' | b'n')),
        "chunked minimizer extraction requires an unambiguous sequence"
    );

    let packed_sequence: PackedNSeqVec = PackedNSeqVec::from_ascii(sequence);
    let packed_sequence_slice = packed_sequence.as_slice();
    let sequence_slice = packed_sequence_slice.seq;
    let hasher: NtHasher<true> = NtHasher::<true>::new_with_seed(kmer_size, minimizer_hash_seed);
    let mut minimizer_positions: Vec<u32> = Vec::new();
    let mut super_kmer_starts: Vec<u32> = Vec::new();
    let _ = canonical_minimizers(kmer_size, window_size)
        .hasher(&hasher)
        .super_kmers(&mut super_kmer_starts)
        .run(sequence_slice, &mut minimizer_positions);
    debug_assert_eq!(minimizer_positions.len(), super_kmer_starts.len());

    minimizer_positions
        .into_iter()
        .zip(super_kmer_starts)
        .filter_map(|(position, super_kmer_start)| {
            let pos: usize = position as usize;
            let forward_kmer = sequence_slice.read_kmer(kmer_size, pos);
            let reverse_complement_kmer = sequence_slice.read_revcomp_kmer(kmer_size, pos);
            if forward_kmer == reverse_complement_kmer {
                return None;
            }

            let canonical_kmer: u64 = forward_kmer.min(reverse_complement_kmer);
            let key: MinimizerKey = MinimizerKey::try_from(canonical_kmer)
                .expect("u32 2-bit minimizer keys require k <= 16");
            Some((key, position, super_kmer_start))
        })
        .collect()
}

fn is_all_n_sequence(sequence: &[u8]) -> bool {
    !sequence.is_empty() && sequence.iter().all(|base| matches!(base, b'N' | b'n'))
}

fn is_ambiguous_base(base: u8) -> bool {
    matches!(base, b'N' | b'n')
}

pub(crate) fn is_no_usable_fragments_error(error: &io::Error) -> bool {
    error.kind() == io::ErrorKind::InvalidData
        && error.to_string() == "ERROR: Input has no usable fragments"
}

/// Return the minimizer hashes used for candidate discovery.
pub(crate) fn select_seed_minimizers(
    minimizers: &[MinimizerKey],
    minmer_count: Option<usize>,
) -> Vec<MinimizerKey> {
    match minmer_count {
        Some(count) => minimizers.iter().take(count).copied().collect(),
        None => minimizers.to_vec(),
    }
}

pub(crate) fn query_fragment_sketch(
    fragment_sequence: &[u8],
    kmer_size: usize,
    window_size: usize,
    minimizer_hash_seed: u32,
    minmer_count: Option<usize>,
) -> QueryFragmentSketch {
    let observation: MinimizerObservation = canonical_minimizer_observation(
        fragment_sequence,
        kmer_size,
        window_size,
        minimizer_hash_seed,
    );
    let mut minimizers: Vec<MinimizerKey> = observation
        .minimizers_with_positions
        .into_iter()
        .map(|(minimizer, _position)| minimizer)
        .collect();
    minimizers.sort_unstable();
    minimizers.dedup();
    let seed_minimizers: Vec<MinimizerKey> = select_seed_minimizers(&minimizers, minmer_count);

    QueryFragmentSketch {
        minimizers,
        seed_minimizers,
    }
}

pub(crate) fn split_sequence_ranges(
    sequence: &[u8],
    split_n_run: usize,
) -> Vec<std::ops::Range<usize>> {
    if sequence.is_empty() {
        return Vec::new();
    }

    if split_n_run == 0 {
        let whole_sequence = 0..sequence.len();
        return vec![whole_sequence];
    }

    let mut ranges: Vec<std::ops::Range<usize>> = Vec::new();
    let mut segment_start: usize = 0usize;
    let mut position: usize = 0usize;

    while position < sequence.len() {
        if !is_ambiguous_base(sequence[position]) {
            position += 1;
            continue;
        }

        let run_start: usize = position;
        while position < sequence.len() && is_ambiguous_base(sequence[position]) {
            position += 1;
        }
        let run_end: usize = position;

        if run_end - run_start >= split_n_run {
            if segment_start < run_start {
                ranges.push(segment_start..run_start);
            }
            segment_start = run_end;
        }
    }

    if segment_start < sequence.len() {
        ranges.push(segment_start..sequence.len());
    }

    ranges
}

/// Return query fragment ranges for either FastANI-compatible chunks or adaptive overlap mode.
pub(crate) fn query_fragment_ranges(
    sequence_len: usize,
    fragment_length: usize,
    fragment_stride: usize,
    min_fragment_length: usize,
) -> Vec<std::ops::Range<usize>> {
    if sequence_len < min_fragment_length {
        return Vec::new();
    }

    let retain_tail_fragment: bool = min_fragment_length < fragment_length;
    let fastani_compatible: bool = fragment_stride == fragment_length && !retain_tail_fragment;

    if fastani_compatible {
        return (0..sequence_len)
            .step_by(fragment_length)
            .take_while(|start| start + fragment_length <= sequence_len)
            .map(|start| start..start + fragment_length)
            .collect();
    }

    if sequence_len <= fragment_length {
        let whole_sequence = 0..sequence_len;
        return vec![whole_sequence];
    }

    let mut ranges: Vec<std::ops::Range<usize>> = Vec::new();
    let mut start: usize = 0usize;

    while start + fragment_length <= sequence_len {
        ranges.push(start..start + fragment_length);
        start += fragment_stride;
    }

    let tail_len: usize = sequence_len.saturating_sub(start);
    if tail_len >= min_fragment_length {
        ranges.push(start..sequence_len);
        return ranges;
    }

    if retain_tail_fragment {
        return ranges;
    }

    let tail_start: usize = sequence_len - fragment_length;
    if ranges
        .last()
        .is_none_or(|last_range| last_range.start != tail_start)
    {
        ranges.push(tail_start..sequence_len);
    }

    ranges
}

pub(crate) fn mapped_length_from_fragment_ranges(
    sequence_len: usize,
    fragment_length: u32,
    min_fragment_length: u32,
) -> u64 {
    query_fragment_ranges(
        sequence_len,
        fragment_length as usize,
        fragment_length as usize,
        min_fragment_length as usize,
    )
    .into_iter()
    .map(|range| range.len() as u64)
    .sum()
}

pub(crate) fn fastani_compatible_fragment_mode(
    fragment_length: u32,
    fragment_stride: u32,
    min_fragment_length: u32,
) -> bool {
    fragment_stride == fragment_length && min_fragment_length == fragment_length
}

#[cfg(test)]
mod tests {
    use crate::ani::{
        constants::{
            MinimizerKey, DEFAULT_FRAGMENT_LENGTH, DEFAULT_FRAGMENT_STRIDE,
            DEFAULT_MINIMIZER_HASH_SEED, DEFAULT_MIN_FRAGMENT_LENGTH,
        },
        io_util::FastaInput,
        minimizer::{
            canonical_minimizer_observation, canonical_minimizers_with_positions,
            expected_minimizer_window_count, fastani_compatible_fragment_mode,
            mapped_length_from_fragment_ranges, query_fragment_ranges, select_seed_minimizers,
            split_sequence_ranges, usable_minimizer_window_count, MinimizerObservation,
        },
        sketch::partition::{
            estimate_reference_minimizer_windows, estimate_selected_minimizers_from_windows,
        },
        test_support::repeated_acgt,
    };
    use std::{env, fs, io, path::PathBuf, time::Instant};

    #[test]
    fn reference_minimizer_window_estimate_respects_split_n() -> io::Result<()> {
        let path: PathBuf = env::temp_dir().join(format!(
            "fasterani_estimate_windows_{}_{}.fa",
            std::process::id(),
            Instant::now().elapsed().as_nanos()
        ));
        fs::write(&path, b">seq\nAAAAAANNNNAAAAAA\n")?;

        let reference = FastaInput::from_path(path.to_str().expect("utf8 temp path").to_string());
        let unsplit_estimate: usize = estimate_reference_minimizer_windows(&reference, 3, 3, 0)?;
        let split_estimate: usize = estimate_reference_minimizer_windows(&reference, 3, 3, 4)?;
        fs::remove_file(&path)?;

        assert_eq!(unsplit_estimate, 12);
        assert_eq!(split_estimate, 4);
        assert_eq!(estimate_selected_minimizers_from_windows(12, 3), 6);
        assert_eq!(estimate_selected_minimizers_from_windows(4, 3), 2);

        Ok(())
    }

    #[test]
    fn minimizers_skip_kmers_with_ambiguous_bases() {
        let sequence: &[u8] =
            b"AGCTTAGGCTAACCGTATGCCGATTAACGNNNNNNNNNNGCTAGTCCATGATCGTACCGTTAAGGCTA";
        let kmer_size: usize = 10usize;
        let window_size: usize = 16usize;

        let minimizers: Vec<(MinimizerKey, u32)> = canonical_minimizers_with_positions(
            sequence,
            kmer_size,
            window_size,
            DEFAULT_MINIMIZER_HASH_SEED,
        );

        assert!(!minimizers.is_empty());
        for (_hash, position) in minimizers {
            let position: usize = position as usize;
            let kmer: &[u8] = &sequence[position..position + kmer_size];
            assert!(!kmer.iter().any(|base| matches!(base, b'N' | b'n')));
        }
    }

    #[test]
    fn all_ambiguous_sequence_has_no_minimizers() {
        let sequence: &[u8] = b"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
        let kmer_size: usize = 10usize;
        let window_size: usize = 16usize;

        let minimizers: Vec<(MinimizerKey, u32)> = canonical_minimizers_with_positions(
            sequence,
            kmer_size,
            window_size,
            DEFAULT_MINIMIZER_HASH_SEED,
        );

        assert!(minimizers.is_empty());
    }

    #[test]
    fn expected_minimizer_windows_match_fragment_geometry() {
        assert_eq!(expected_minimizer_window_count(3000, 16, 24), 2962);
        assert_eq!(expected_minimizer_window_count(38, 16, 24), 0);
        assert_eq!(expected_minimizer_window_count(39, 16, 24), 1);
    }

    #[test]
    fn default_fragment_settings_are_fastani_compatible() {
        assert_eq!(DEFAULT_FRAGMENT_STRIDE, DEFAULT_FRAGMENT_LENGTH);
        assert_eq!(DEFAULT_MIN_FRAGMENT_LENGTH, DEFAULT_FRAGMENT_LENGTH);
        assert!(fastani_compatible_fragment_mode(
            DEFAULT_FRAGMENT_LENGTH,
            DEFAULT_FRAGMENT_STRIDE,
            DEFAULT_MIN_FRAGMENT_LENGTH
        ));
    }

    #[test]
    fn clean_sequence_has_no_missing_minimizer_windows() {
        let sequence: Vec<u8> = repeated_acgt(120);

        assert_eq!(
            usable_minimizer_window_count(&sequence, 5, 5),
            expected_minimizer_window_count(sequence.len(), 5, 5)
        );
    }

    #[test]
    fn all_ambiguous_observation_keeps_expected_window_count() {
        let sequence: &[u8] = b"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
        let observation: MinimizerObservation =
            canonical_minimizer_observation(sequence, 10, 16, DEFAULT_MINIMIZER_HASH_SEED);

        assert_eq!(usable_minimizer_window_count(sequence, 10, 16), 0);
        assert_eq!(expected_minimizer_window_count(sequence.len(), 10, 16), 16);
        assert!(observation.minimizers_with_positions.is_empty());
    }

    #[test]
    fn minmer_seeds_use_smallest_hashes() {
        let minimizers: Vec<MinimizerKey> = vec![3, 5, 8, 13, 21];

        let seeds: Vec<MinimizerKey> = select_seed_minimizers(&minimizers, Some(3));

        assert_eq!(seeds, vec![3, 5, 8]);
    }

    #[test]
    fn default_fragment_ranges_match_fastani_chunks() {
        let ranges: Vec<std::ops::Range<usize>> = query_fragment_ranges(7500, 3000, 3000, 3000);

        assert_eq!(ranges, vec![0..3000, 3000..6000]);
    }

    #[test]
    fn default_fragment_ranges_discard_short_terminal_remainder() {
        let ranges: Vec<std::ops::Range<usize>> = query_fragment_ranges(8999, 3000, 3000, 3000);

        assert_eq!(ranges, vec![0..3000, 3000..6000]);
    }

    #[test]
    fn adaptive_fragment_ranges_overlap_and_cover_tail() {
        let ranges: Vec<std::ops::Range<usize>> = query_fragment_ranges(7500, 3000, 1000, 1000);

        assert_eq!(
            ranges,
            vec![
                0..3000,
                1000..4000,
                2000..5000,
                3000..6000,
                4000..7000,
                5000..7500
            ]
        );
    }

    #[test]
    fn adaptive_fragment_ranges_keep_short_usable_contigs() {
        let ranges: Vec<std::ops::Range<usize>> = query_fragment_ranges(1500, 3000, 1000, 1000);

        assert_eq!(ranges, vec![0..1500]);
    }

    #[test]
    fn adaptive_fragment_ranges_keep_terminal_remainder() {
        let ranges: Vec<std::ops::Range<usize>> = query_fragment_ranges(7500, 3000, 3000, 1000);

        assert_eq!(ranges, vec![0..3000, 3000..6000, 6000..7500]);
    }

    #[test]
    fn mapped_length_keeps_terminal_remainder_when_allowed() {
        assert_eq!(
            mapped_length_from_fragment_ranges(7500, DEFAULT_FRAGMENT_LENGTH, 3000),
            6000
        );
        assert_eq!(
            mapped_length_from_fragment_ranges(7500, DEFAULT_FRAGMENT_LENGTH, 1000),
            7500
        );
        assert_eq!(
            mapped_length_from_fragment_ranges(3500, DEFAULT_FRAGMENT_LENGTH, 1000),
            3000
        );
        assert_eq!(
            mapped_length_from_fragment_ranges(4000, DEFAULT_FRAGMENT_LENGTH, 1000),
            4000
        );
    }

    #[test]
    fn split_n_zero_keeps_whole_sequence() {
        let ranges: Vec<std::ops::Range<usize>> = split_sequence_ranges(b"ACGTNNNNACGT", 0);

        assert_eq!(ranges, vec![0..12]);
    }

    #[test]
    fn split_n_keeps_short_ambiguous_runs() {
        let ranges: Vec<std::ops::Range<usize>> = split_sequence_ranges(b"ACGTNNACGT", 3);

        assert_eq!(ranges, vec![0..10]);
    }

    #[test]
    fn split_n_breaks_long_ambiguous_runs() {
        let ranges: Vec<std::ops::Range<usize>> = split_sequence_ranges(b"NNNACGTNNNNACGTNNN", 3);

        assert_eq!(ranges, vec![3..7, 11..15]);
    }
}
