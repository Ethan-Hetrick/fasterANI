//! Database-global minimizer frequencies for representation-independent filtering.

use std::{
    cmp::Reverse,
    collections::{BinaryHeap, HashMap},
    fs, io,
    io::{BufReader, BufWriter, Read, Write},
    path::{Path, PathBuf},
    sync::Arc,
    time::Instant,
};

use crate::ani::{
    io_util::ScratchFile,
    mmap::MmapFile,
    model::reference::{
        ReferenceIndex, ReferenceSketch, ShardManifest, ShardManifestEntry, SketchParams,
    },
    runtime::{emit_progress, RuntimeOptions},
    sketch::serialize::{
        global_frequency_entry_path, global_frequency_filename, global_frequency_path, SketchOutput,
    },
};

const GLOBAL_FREQUENCY_MAGIC: &[u8; 8] = b"FANIFRQ1";
const GLOBAL_FREQUENCY_FORMAT_VERSION: u32 = 2;
const FREQUENCY_RECORD_BYTES: usize = 12;
const HISTOGRAM_RECORD_BYTES: usize = 16;
const MERGE_FAN_IN: usize = 64;

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct FrequencyRecord {
    key: u32,
    count: u64,
}

struct FrequencyRun {
    scratch: ScratchFile,
    record_count: u64,
}

struct RunReader {
    reader: BufReader<fs::File>,
    remaining: u64,
}

impl RunReader {
    fn open(run: &FrequencyRun) -> io::Result<Self> {
        let expected_bytes: u64 = run
            .record_count
            .checked_mul(FREQUENCY_RECORD_BYTES as u64)
            .ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "frequency run size overflow")
            })?;
        let actual_bytes: u64 = fs::metadata(&run.scratch.path)?.len();
        if actual_bytes != expected_bytes {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "global-frequency scratch run {} has {actual_bytes} bytes, expected {expected_bytes}",
                    run.scratch.path.display()
                ),
            ));
        }
        Ok(Self {
            reader: BufReader::new(fs::File::open(&run.scratch.path)?),
            remaining: run.record_count,
        })
    }

    fn next_record(&mut self) -> io::Result<Option<FrequencyRecord>> {
        if self.remaining == 0 {
            return Ok(None);
        }
        let record: FrequencyRecord = read_frequency_record(&mut self.reader)?;
        self.remaining -= 1;
        Ok(Some(record))
    }
}

pub(crate) struct GlobalFrequencyArtifactStats {
    pub(crate) filename: String,
    pub(crate) file_bytes: u64,
    pub(crate) unique_minimizers: usize,
}

/// Read-only database-global key frequencies and their exact frequency histogram.
pub(crate) struct GlobalFrequencyIndex {
    mmap: Arc<MmapFile>,
    records_offset: usize,
    record_count: usize,
    histogram_offset: usize,
    histogram_count: usize,
}

impl GlobalFrequencyIndex {
    pub(crate) fn load(prefix: &Path, manifest: &ShardManifest) -> io::Result<Self> {
        if manifest.global_frequency_filename.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "sharded sketch manifest has no database-global frequency artifact; rebuild the sketch database",
            ));
        }
        let expected_filename: String = global_frequency_filename(prefix, &manifest.generation_id);
        if manifest.global_frequency_filename != expected_filename {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "global-frequency artifact filename mismatch: manifest={} expected={expected_filename}; rebuild the sketch database",
                    manifest.global_frequency_filename
                ),
            ));
        }
        let path: PathBuf =
            global_frequency_entry_path(prefix, &manifest.global_frequency_filename);
        let actual_bytes: u64 = fs::metadata(&path)
            .map_err(|err| {
                io::Error::new(
                    err.kind(),
                    format!(
                        "sharded sketch manifest references unreadable global-frequency artifact {}: {err}",
                        path.display()
                    ),
                )
            })?
            .len();
        if actual_bytes != manifest.global_frequency_file_bytes {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "global-frequency artifact size mismatch for {}: manifest={} actual={}; rebuild the sketch database",
                    path.display(),
                    manifest.global_frequency_file_bytes,
                    actual_bytes
                ),
            ));
        }

        let mmap: Arc<MmapFile> = Arc::new(MmapFile::open(&path)?);
        let bytes: &[u8] = mmap.as_slice();
        let fixed_header_bytes: usize = GLOBAL_FREQUENCY_MAGIC.len() + 4 + 8 + 8;
        if bytes.len() < fixed_header_bytes || &bytes[..8] != GLOBAL_FREQUENCY_MAGIC {
            return Err(invalid_artifact(&path, "invalid magic header"));
        }
        let format_version: u32 = read_u32_at(bytes, 8, &path)?;
        if format_version != GLOBAL_FREQUENCY_FORMAT_VERSION {
            return Err(invalid_artifact(
                &path,
                &format!(
                    "format version {format_version} is not supported; rebuild the sketch database"
                ),
            ));
        }
        let record_count: usize = usize::try_from(read_u64_at(bytes, 12, &path)?)
            .map_err(|_| invalid_artifact(&path, "record count does not fit this platform"))?;
        let histogram_count: usize = usize::try_from(read_u64_at(bytes, 20, &path)?)
            .map_err(|_| invalid_artifact(&path, "histogram count does not fit this platform"))?;
        let records_offset: usize = fixed_header_bytes;
        let histogram_offset: usize = records_offset
            .checked_add(
                record_count
                    .checked_mul(FREQUENCY_RECORD_BYTES)
                    .ok_or_else(|| invalid_artifact(&path, "frequency record size overflow"))?,
            )
            .ok_or_else(|| invalid_artifact(&path, "frequency section offset overflow"))?;
        let expected_file_bytes: usize = histogram_offset
            .checked_add(
                histogram_count
                    .checked_mul(HISTOGRAM_RECORD_BYTES)
                    .ok_or_else(|| invalid_artifact(&path, "histogram size overflow"))?,
            )
            .ok_or_else(|| invalid_artifact(&path, "artifact size overflow"))?;
        if expected_file_bytes != bytes.len() {
            return Err(invalid_artifact(
                &path,
                &format!(
                    "array sections require {expected_file_bytes} bytes but the file has {}",
                    bytes.len()
                ),
            ));
        }
        if record_count != manifest.total_unique_minimizers {
            return Err(invalid_artifact(
                &path,
                &format!(
                    "record count {record_count} does not match manifest total {}",
                    manifest.total_unique_minimizers
                ),
            ));
        }

        let index = Self {
            mmap,
            records_offset,
            record_count,
            histogram_offset,
            histogram_count,
        };
        index.validate_sorted_contents(&path)?;
        Ok(index)
    }

    pub(crate) fn get(&self, key: u32) -> Option<usize> {
        let mut left: usize = 0;
        let mut right: usize = self.record_count;
        while left < right {
            let middle: usize = left + (right - left) / 2;
            let (middle_key, count): (u32, u64) = self.record(middle);
            match middle_key.cmp(&key) {
                std::cmp::Ordering::Less => left = middle + 1,
                std::cmp::Ordering::Greater => right = middle,
                std::cmp::Ordering::Equal => {
                    return Some(usize::try_from(count).unwrap_or(usize::MAX));
                }
            }
        }
        None
    }

    pub(crate) fn frequency_threshold(&self, percent: f64) -> usize {
        if percent <= 0.0 {
            return usize::MAX;
        }
        let minimizers_to_ignore: usize = (self.record_count as f64 * percent / 100.0) as usize;
        if minimizers_to_ignore == 0 {
            return usize::MAX;
        }

        let mut sum: usize = 0;
        let mut threshold: usize = usize::MAX;
        for index in 0..self.histogram_count {
            let (frequency, key_count): (u64, u64) = self.histogram_record(index);
            let frequency: usize = usize::try_from(frequency).unwrap_or(usize::MAX);
            let key_count: usize = usize::try_from(key_count).unwrap_or(usize::MAX);
            sum = sum.saturating_add(key_count);
            match sum.cmp(&minimizers_to_ignore) {
                std::cmp::Ordering::Less => threshold = frequency,
                std::cmp::Ordering::Equal => {
                    threshold = frequency;
                    break;
                }
                std::cmp::Ordering::Greater => break,
            }
        }
        threshold
    }

    fn record(&self, index: usize) -> (u32, u64) {
        let start: usize = self.records_offset + index * FREQUENCY_RECORD_BYTES;
        let bytes: &[u8] = self.mmap.as_slice();
        (
            u32::from_le_bytes(
                bytes[start..start + 4]
                    .try_into()
                    .expect("validated record"),
            ),
            u64::from_le_bytes(
                bytes[start + 4..start + 12]
                    .try_into()
                    .expect("validated record"),
            ),
        )
    }

    fn histogram_record(&self, index: usize) -> (u64, u64) {
        let start: usize = self.histogram_offset + index * HISTOGRAM_RECORD_BYTES;
        let bytes: &[u8] = self.mmap.as_slice();
        (
            u64::from_le_bytes(
                bytes[start..start + 8]
                    .try_into()
                    .expect("validated histogram"),
            ),
            u64::from_le_bytes(
                bytes[start + 8..start + 16]
                    .try_into()
                    .expect("validated histogram"),
            ),
        )
    }

    fn validate_sorted_contents(&self, path: &Path) -> io::Result<()> {
        let mut previous_key: Option<u32> = None;
        for index in 0..self.record_count {
            let (key, count): (u32, u64) = self.record(index);
            if count == 0 || previous_key.is_some_and(|previous| previous >= key) {
                return Err(invalid_artifact(
                    path,
                    "frequency records are not strictly key-sorted with positive counts",
                ));
            }
            previous_key = Some(key);
        }

        let mut previous_frequency: Option<u64> = None;
        let mut histogram_keys: u64 = 0;
        for index in 0..self.histogram_count {
            let (frequency, key_count): (u64, u64) = self.histogram_record(index);
            if frequency == 0
                || key_count == 0
                || previous_frequency.is_some_and(|previous| previous <= frequency)
            {
                return Err(invalid_artifact(
                    path,
                    "frequency histogram is not strictly descending with positive counts",
                ));
            }
            histogram_keys = histogram_keys
                .checked_add(key_count)
                .ok_or_else(|| invalid_artifact(path, "frequency histogram key count overflow"))?;
            previous_frequency = Some(frequency);
        }
        if histogram_keys != self.record_count as u64 {
            return Err(invalid_artifact(
                path,
                "frequency histogram does not cover every global minimizer key",
            ));
        }
        Ok(())
    }
}

pub(crate) fn build_global_frequency_artifact(
    prefix: &Path,
    generation_id: &str,
    shards: &[ShardManifestEntry],
    params: SketchParams,
    tmp_dir: Option<&Path>,
    runtime_options: RuntimeOptions,
) -> io::Result<GlobalFrequencyArtifactStats> {
    let build_start: Instant = Instant::now();
    if runtime_options.progress_enabled {
        emit_progress(
            "database_build",
            &format!(
                "event=global_frequency_start\tgeneration_id={generation_id}\tshards={}",
                shards.len()
            ),
            build_start,
        );
    }
    let mut runs: Vec<FrequencyRun> = Vec::with_capacity(shards.len());
    for (shard_offset, shard) in shards.iter().enumerate() {
        let shard_runtime_options: RuntimeOptions =
            runtime_options.with_build_progress(generation_id, shard.shard_index)?;
        let sketch: ReferenceSketch = ReferenceSketch::load(
            &crate::ani::sketch::serialize::shard_entry_path(prefix, shard),
            params,
            false,
            shard_runtime_options,
        )?;
        let ReferenceIndex::Mphf(index) = &sketch.index else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "persisted sketch shard did not load as an MPHF index",
            ));
        };
        let mut records: Vec<FrequencyRecord> = index
            .slot_keys()
            .iter()
            .copied()
            .zip(index.hit_counts().iter().copied())
            .map(|(key, count)| FrequencyRecord {
                key,
                count: u64::from(count),
            })
            .collect();
        records.sort_unstable_by_key(|record| record.key);
        runs.push(write_run(&records, tmp_dir, "global-frequency-shard")?);
        if runtime_options.progress_enabled {
            emit_progress(
                "database_build",
                &format!(
                    "event=global_frequency_shard\tgeneration_id={generation_id}\tshard={}\tshards_done={}\tshards_total={}\tlocal_unique_minimizers={}",
                    shard.shard_index,
                    shard_offset + 1,
                    shards.len(),
                    records.len()
                ),
                build_start,
            );
        }
    }

    let merged_run: Option<FrequencyRun> = merge_all_runs(runs, tmp_dir)?;
    let final_path: PathBuf = global_frequency_path(prefix, generation_id);
    let (file_bytes, unique_minimizers): (u64, usize) =
        write_frequency_artifact(&final_path, merged_run.as_ref())?;
    let stats = GlobalFrequencyArtifactStats {
        filename: global_frequency_filename(prefix, generation_id),
        file_bytes,
        unique_minimizers,
    };
    if runtime_options.progress_enabled {
        emit_progress(
            "database_build",
            &format!(
                "event=global_frequency_complete\tgeneration_id={generation_id}\tunique_minimizers={}\tfile_bytes={}\tfilename={}",
                stats.unique_minimizers, stats.file_bytes, stats.filename
            ),
            build_start,
        );
    }
    Ok(stats)
}

fn write_run(
    records: &[FrequencyRecord],
    tmp_dir: Option<&Path>,
    purpose: &str,
) -> io::Result<FrequencyRun> {
    let (scratch, file): (ScratchFile, fs::File) = ScratchFile::create(tmp_dir, purpose)?;
    let mut writer: BufWriter<fs::File> = BufWriter::new(file);
    for &record in records {
        write_frequency_record(&mut writer, record)?;
    }
    writer.flush()?;
    Ok(FrequencyRun {
        scratch,
        record_count: records.len() as u64,
    })
}

fn merge_all_runs(
    mut runs: Vec<FrequencyRun>,
    tmp_dir: Option<&Path>,
) -> io::Result<Option<FrequencyRun>> {
    while runs.len() > 1 {
        let mut next_runs: Vec<FrequencyRun> =
            Vec::with_capacity(runs.len().div_ceil(MERGE_FAN_IN));
        let mut remaining = runs.into_iter();
        loop {
            let group: Vec<FrequencyRun> = remaining.by_ref().take(MERGE_FAN_IN).collect();
            if group.is_empty() {
                break;
            }
            if group.len() == 1 {
                next_runs.push(group.into_iter().next().expect("one run"));
            } else {
                next_runs.push(merge_run_group(&group, tmp_dir)?);
            }
        }
        runs = next_runs;
    }
    Ok(runs.pop())
}

fn merge_run_group(runs: &[FrequencyRun], tmp_dir: Option<&Path>) -> io::Result<FrequencyRun> {
    let mut readers: Vec<RunReader> = runs
        .iter()
        .map(RunReader::open)
        .collect::<io::Result<Vec<_>>>()?;
    let mut heap: BinaryHeap<Reverse<(u32, usize, u64)>> = BinaryHeap::new();
    for (run_index, reader) in readers.iter_mut().enumerate() {
        if let Some(record) = reader.next_record()? {
            heap.push(Reverse((record.key, run_index, record.count)));
        }
    }

    let (scratch, file): (ScratchFile, fs::File) =
        ScratchFile::create(tmp_dir, "global-frequency-merge")?;
    let mut writer: BufWriter<fs::File> = BufWriter::new(file);
    let mut output_count: u64 = 0;
    let mut pending: Option<FrequencyRecord> = None;
    while let Some(Reverse((key, run_index, count))) = heap.pop() {
        if let Some(next_record) = readers[run_index].next_record()? {
            heap.push(Reverse((next_record.key, run_index, next_record.count)));
        }
        match pending.as_mut() {
            Some(record) if record.key == key => {
                record.count = record.count.checked_add(count).ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        "global minimizer count overflow",
                    )
                })?;
            }
            Some(record) => {
                write_frequency_record(&mut writer, *record)?;
                output_count += 1;
                *record = FrequencyRecord { key, count };
            }
            None => pending = Some(FrequencyRecord { key, count }),
        }
    }
    if let Some(record) = pending {
        write_frequency_record(&mut writer, record)?;
        output_count += 1;
    }
    writer.flush()?;
    Ok(FrequencyRun {
        scratch,
        record_count: output_count,
    })
}

fn write_frequency_artifact(
    path: &Path,
    merged_run: Option<&FrequencyRun>,
) -> io::Result<(u64, usize)> {
    let record_count: u64 = merged_run.map_or(0, |run| run.record_count);
    let mut histogram: HashMap<u64, u64> = HashMap::new();
    if let Some(run) = merged_run {
        let mut reader: RunReader = RunReader::open(run)?;
        while let Some(record) = reader.next_record()? {
            *histogram.entry(record.count).or_default() += 1;
        }
    }
    let mut histogram: Vec<(u64, u64)> = histogram.into_iter().collect();
    histogram.sort_unstable_by_key(|&(frequency, _)| Reverse(frequency));
    let mut output: SketchOutput = SketchOutput::create(path)?;
    let writer: &mut BufWriter<fs::File> = output.writer_mut()?;
    writer.write_all(GLOBAL_FREQUENCY_MAGIC)?;
    writer.write_all(&GLOBAL_FREQUENCY_FORMAT_VERSION.to_le_bytes())?;
    writer.write_all(&record_count.to_le_bytes())?;
    writer.write_all(&(histogram.len() as u64).to_le_bytes())?;
    if let Some(run) = merged_run {
        let mut reader: BufReader<fs::File> = BufReader::new(fs::File::open(&run.scratch.path)?);
        io::copy(&mut reader, writer)?;
    }
    for (frequency, key_count) in histogram {
        writer.write_all(&frequency.to_le_bytes())?;
        writer.write_all(&key_count.to_le_bytes())?;
    }
    let file_bytes: u64 = output.finish()?;
    let unique_minimizers: usize = usize::try_from(record_count).map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            "global unique-minimizer count does not fit this platform",
        )
    })?;
    Ok((file_bytes, unique_minimizers))
}

fn write_frequency_record(writer: &mut impl Write, record: FrequencyRecord) -> io::Result<()> {
    writer.write_all(&record.key.to_le_bytes())?;
    writer.write_all(&record.count.to_le_bytes())
}

fn read_frequency_record(reader: &mut impl Read) -> io::Result<FrequencyRecord> {
    let mut key: [u8; 4] = [0; 4];
    let mut count: [u8; 8] = [0; 8];
    reader.read_exact(&mut key)?;
    reader.read_exact(&mut count)?;
    Ok(FrequencyRecord {
        key: u32::from_le_bytes(key),
        count: u64::from_le_bytes(count),
    })
}

fn read_u32_at(bytes: &[u8], offset: usize, path: &Path) -> io::Result<u32> {
    let value: [u8; 4] = bytes
        .get(offset..offset + 4)
        .ok_or_else(|| invalid_artifact(path, "truncated header"))?
        .try_into()
        .expect("slice length checked");
    Ok(u32::from_le_bytes(value))
}

fn read_u64_at(bytes: &[u8], offset: usize, path: &Path) -> io::Result<u64> {
    let value: [u8; 8] = bytes
        .get(offset..offset + 8)
        .ok_or_else(|| invalid_artifact(path, "truncated header"))?
        .try_into()
        .expect("slice length checked");
    Ok(u64::from_le_bytes(value))
}

fn invalid_artifact(path: &Path, detail: &str) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidData,
        format!(
            "invalid global-frequency artifact {}: {detail}",
            path.display()
        ),
    )
}

#[cfg(test)]
mod tests {
    use super::{
        merge_all_runs, write_frequency_artifact, write_run, FrequencyRecord, GlobalFrequencyIndex,
    };
    use crate::ani::{
        constants::ReferenceHitMap,
        model::reference::{ReferenceContigs, ReferenceIndex, ReferenceSketch, SeedHit},
        sketch::serialize::global_frequency_filename,
    };
    use std::{env, fs, io, path::PathBuf, sync::Arc, time::SystemTime};

    #[test]
    fn distributed_key_uses_global_count_and_monolithic_threshold() -> io::Result<()> {
        let first = write_run(
            &[
                FrequencyRecord { key: 7, count: 2 },
                FrequencyRecord { key: 11, count: 1 },
            ],
            None,
            "global-frequency-test-first",
        )?;
        let second = write_run(
            &[
                FrequencyRecord { key: 7, count: 3 },
                FrequencyRecord { key: 13, count: 4 },
            ],
            None,
            "global-frequency-test-second",
        )?;
        let merged = merge_all_runs(vec![first, second], None)?.expect("merged run");
        let unique: u128 = SystemTime::now()
            .duration_since(SystemTime::UNIX_EPOCH)
            .map_err(io::Error::other)?
            .as_nanos();
        let prefix: PathBuf = env::temp_dir().join(format!(
            "fasterani-global-frequency-test-{}-{unique}",
            std::process::id()
        ));
        let generation = "test-generation-one";
        let second_generation = "test-generation-two";
        let artifact_path: PathBuf = super::global_frequency_path(&prefix, generation);
        let second_artifact_path: PathBuf =
            super::global_frequency_path(&prefix, second_generation);
        let (file_bytes, unique_minimizers) =
            write_frequency_artifact(&artifact_path, Some(&merged))?;
        let (second_file_bytes, second_unique_minimizers) =
            write_frequency_artifact(&second_artifact_path, Some(&merged))?;
        assert_eq!(file_bytes, second_file_bytes);
        assert_eq!(unique_minimizers, second_unique_minimizers);
        assert_eq!(
            fs::read(&artifact_path)?,
            fs::read(&second_artifact_path)?,
            "global-frequency bytes must not depend on the build generation"
        );
        let mut manifest = crate::ani::test_support::sample_shard_manifest();
        manifest.generation_id = generation.to_string();
        manifest.global_frequency_filename = global_frequency_filename(&prefix, generation);
        manifest.global_frequency_file_bytes = file_bytes;
        manifest.total_unique_minimizers = unique_minimizers;
        let global = Arc::new(GlobalFrequencyIndex::load(&prefix, &manifest)?);

        assert_eq!(global.get(7), Some(5));
        assert_eq!(global.get(11), Some(1));
        assert_eq!(global.get(13), Some(4));

        let mut monolithic: ReferenceHitMap = ReferenceHitMap::default();
        monolithic.insert(7, vec![SeedHit::default(); 5]);
        monolithic.insert(11, vec![SeedHit::default(); 1]);
        monolithic.insert(13, vec![SeedHit::default(); 4]);
        let monolithic = ReferenceIndex::Hash(monolithic);
        for percent in [0.0, 10.0, 33.4, 34.0, 66.7, 100.0] {
            assert_eq!(
                global.frequency_threshold(percent),
                monolithic.frequency_threshold(percent),
                "threshold mismatch at {percent}%"
            );
        }

        let mut local_hits: ReferenceHitMap = ReferenceHitMap::default();
        local_hits.insert(
            7,
            vec![
                SeedHit {
                    reference_contig_id: 0,
                    position: 5,
                },
                SeedHit {
                    reference_contig_id: 0,
                    position: 10,
                },
            ],
        );
        let sketch = ReferenceSketch {
            files: Vec::new(),
            contigs: ReferenceContigs::Owned(Vec::new()),
            contig_names: None,
            index: ReferenceIndex::Hash(local_hits),
            global_frequencies: Some(global),
        };
        let mut seed_hits = Vec::new();
        let mut candidates = Vec::new();
        let mut slots = Vec::new();
        let mut ranges = Vec::new();
        sketch.find_candidate_regions(
            &[7],
            100,
            1,
            5,
            &mut seed_hits,
            &mut candidates,
            &mut slots,
            &mut ranges,
            #[cfg(debug_assertions)]
            None,
        );
        assert!(
            seed_hits.is_empty(),
            "global count must filter the local hits"
        );
        assert!(candidates.is_empty());

        manifest.generation_id = second_generation.to_string();
        let filename_error = match GlobalFrequencyIndex::load(&prefix, &manifest) {
            Ok(_) => {
                panic!("a manifest must not reference another generation's artifact filename")
            }
            Err(error) => error,
        };
        assert!(
            filename_error.to_string().contains("filename mismatch"),
            "unexpected filename validation error: {filename_error}"
        );

        fs::remove_file(artifact_path)?;
        fs::remove_file(second_artifact_path)?;
        Ok(())
    }
}

/// Update global counts by subtracting old affected shards and adding their
/// replacements/new shards. Unchanged shard indexes never need to be rescanned.
#[allow(clippy::too_many_arguments)]
pub(crate) fn update_global_frequency_artifact(
    prefix: &Path,
    generation: &str,
    old: &GlobalFrequencyIndex,
    removed: &[ShardManifestEntry],
    added: &[ShardManifestEntry],
    params: SketchParams,
    tmp_dir: Option<&Path>,
    runtime: RuntimeOptions,
) -> io::Result<GlobalFrequencyArtifactStats> {
    fn shard_runs(
        prefix: &Path,
        shards: &[ShardManifestEntry],
        params: SketchParams,
        tmp: Option<&Path>,
        runtime: RuntimeOptions,
    ) -> io::Result<Option<FrequencyRun>> {
        let mut runs = Vec::new();
        for entry in shards {
            let sketch = ReferenceSketch::load(
                &crate::ani::sketch::serialize::shard_entry_path(prefix, entry),
                params,
                false,
                runtime,
            )?;
            let ReferenceIndex::Mphf(index) = &sketch.index else {
                return Err(io::Error::other("expected persisted MPHF index"));
            };
            let mut records: Vec<_> = index
                .slot_keys()
                .iter()
                .zip(index.hit_counts())
                .map(|(&key, &count)| FrequencyRecord {
                    key,
                    count: u64::from(count),
                })
                .collect();
            records.sort_unstable_by_key(|r| r.key);
            runs.push(write_run(&records, tmp, "update-frequency-delta")?);
        }
        merge_all_runs(runs, tmp)
    }
    let negative = shard_runs(prefix, removed, params, tmp_dir, runtime)?;
    let positive = shard_runs(prefix, added, params, tmp_dir, runtime)?;
    let mut minus = negative.as_ref().map(RunReader::open).transpose()?;
    let mut plus = positive.as_ref().map(RunReader::open).transpose()?;
    fn advance(reader: &mut Option<RunReader>) -> io::Result<Option<FrequencyRecord>> {
        reader
            .as_mut()
            .map(RunReader::next_record)
            .transpose()
            .map(Option::flatten)
    }
    let mut neg = advance(&mut minus)?;
    let mut pos = advance(&mut plus)?;
    let mut old_index = 0;
    let (scratch, file) = ScratchFile::create(tmp_dir, "update-global-frequencies")?;
    let mut writer = BufWriter::new(file);
    let mut record_count = 0;
    loop {
        let previous = (old_index < old.record_count).then(|| old.record(old_index));
        let key = previous
            .map(|r| r.0)
            .into_iter()
            .chain(neg.map(|r| r.key))
            .chain(pos.map(|r| r.key))
            .min();
        let Some(key) = key else {
            break;
        };
        let mut count = 0u64;
        if let Some((old_key, old_count)) = previous {
            if old_key == key {
                count = old_count;
                old_index += 1;
            }
        }
        if neg.is_some_and(|r| r.key == key) {
            count = count.checked_sub(neg.unwrap().count).ok_or_else(|| {
                io::Error::other("removed frequencies exceed original database counts")
            })?;
            neg = advance(&mut minus)?;
        }
        if pos.is_some_and(|r| r.key == key) {
            count = count
                .checked_add(pos.unwrap().count)
                .ok_or_else(|| io::Error::other("frequency count overflow"))?;
            pos = advance(&mut plus)?;
        }
        if count != 0 {
            write_frequency_record(&mut writer, FrequencyRecord { key, count })?;
            record_count += 1;
        }
    }
    writer.flush()?;
    drop(writer);
    let run = FrequencyRun {
        scratch,
        record_count,
    };
    let (file_bytes, unique_minimizers) =
        write_frequency_artifact(&global_frequency_path(prefix, generation), Some(&run))?;
    Ok(GlobalFrequencyArtifactStats {
        filename: global_frequency_filename(prefix, generation),
        file_bytes,
        unique_minimizers,
    })
}
