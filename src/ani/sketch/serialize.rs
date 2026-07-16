//! Reading and writing sketch files, shard manifests, and contig sidecars.

use std::{
    fs,
    fs::OpenOptions,
    io,
    io::{BufReader, BufWriter, Read, Write},
    path::{Component, Path, PathBuf},
    sync::atomic::{AtomicU64, Ordering},
    time::{SystemTime, UNIX_EPOCH},
};

use gzp::{deflate::Bgzf, ZBuilder};

use crate::ani::{
    append_path_suffix, compress_file_to_bgzf, gzp_error_to_io, is_gzip_path, read_text_maybe_gzip,
    ContigRecord, FastaInput, ReferenceContigName, ReferenceFile, ScratchFile, ShardManifestEntry,
};

pub(crate) struct SketchOutput {
    pub(crate) final_path: PathBuf,
    pub(crate) write_path: PathBuf,
    pub(crate) publish_path: PathBuf,
    pub(crate) scratch: Option<ScratchFile>,
    pub(crate) writer: Option<BufWriter<fs::File>>,
    pub(crate) bgzip: bool,
    pub(crate) threads: usize,
    published: bool,
}

static NEXT_TEMP_FILE_ID: AtomicU64 = AtomicU64::new(0);

fn create_sibling_temp_file(final_path: &Path) -> io::Result<(PathBuf, fs::File)> {
    if let Some(parent) = final_path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        fs::create_dir_all(parent)?;
    }

    for _ in 0..1024 {
        let id: u64 = NEXT_TEMP_FILE_ID.fetch_add(1, Ordering::Relaxed);
        let suffix: String = format!(".tmp.{}.{id}", std::process::id());
        let path: PathBuf = append_path_suffix(final_path, &suffix);
        match OpenOptions::new().write(true).create_new(true).open(&path) {
            Ok(file) => return Ok((path, file)),
            Err(err) if err.kind() == io::ErrorKind::AlreadyExists => continue,
            Err(err) => return Err(err),
        }
    }

    Err(io::Error::new(
        io::ErrorKind::AlreadyExists,
        format!(
            "failed to allocate a temporary output next to {}",
            final_path.display()
        ),
    ))
}

fn sync_parent(path: &Path) -> io::Result<()> {
    let parent: &Path = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    fs::File::open(parent)?.sync_all()
}

fn publish_temp_file(temp_path: &Path, final_path: &Path) -> io::Result<()> {
    fs::File::open(temp_path)?.sync_all()?;
    fs::rename(temp_path, final_path)?;
    sync_parent(final_path)
}

pub(crate) fn write_bytes_atomically(path: &Path, bytes: &[u8]) -> io::Result<()> {
    let (temp_path, file): (PathBuf, fs::File) = create_sibling_temp_file(path)?;
    let result: io::Result<()> = (|| {
        let mut writer: BufWriter<fs::File> = BufWriter::new(file);
        writer.write_all(bytes)?;
        writer.flush()?;
        drop(writer);
        publish_temp_file(&temp_path, path)
    })();
    if result.is_err() {
        let _ = fs::remove_file(&temp_path);
    }
    result
}

impl SketchOutput {
    pub(crate) fn create(
        final_path: &Path,
        tmp_dir: Option<&Path>,
        bgzip: bool,
        threads: usize,
    ) -> io::Result<Self> {
        if let Some(parent) = final_path
            .parent()
            .filter(|parent| !parent.as_os_str().is_empty())
        {
            fs::create_dir_all(parent)?;
        }

        let (publish_path, publish_file): (PathBuf, fs::File) =
            create_sibling_temp_file(final_path)?;
        if bgzip {
            let (scratch, file): (ScratchFile, fs::File) =
                ScratchFile::create(tmp_dir, "uncompressed-sketch")?;
            let write_path: PathBuf = scratch.path.clone();
            drop(publish_file);
            Ok(Self {
                final_path: final_path.to_path_buf(),
                write_path,
                publish_path,
                scratch: Some(scratch),
                writer: Some(BufWriter::new(file)),
                bgzip,
                threads,
                published: false,
            })
        } else {
            Ok(Self {
                final_path: final_path.to_path_buf(),
                write_path: publish_path.clone(),
                publish_path,
                scratch: None,
                writer: Some(BufWriter::new(publish_file)),
                bgzip,
                threads,
                published: false,
            })
        }
    }

    pub(crate) fn writer_mut(&mut self) -> io::Result<&mut BufWriter<fs::File>> {
        self.writer
            .as_mut()
            .ok_or_else(|| io::Error::other("reference sketch output writer is already closed"))
    }

    pub(crate) fn finish(mut self) -> io::Result<u64> {
        if let Some(mut writer) = self.writer.take() {
            writer.flush()?;
        }

        if self.bgzip {
            compress_file_to_bgzf(&self.write_path, &self.publish_path, self.threads)?;
        }

        publish_temp_file(&self.publish_path, &self.final_path)?;
        self.published = true;
        let output_len: u64 = fs::metadata(&self.final_path)?.len();
        drop(self.scratch.take());
        Ok(output_len)
    }
}

impl Drop for SketchOutput {
    fn drop(&mut self) {
        if !self.published {
            let _ = fs::remove_file(&self.publish_path);
        }
    }
}

fn content_addressed_contig_sidecar_path(sketch_path: &Path, content_id: &str) -> PathBuf {
    let filename: String = if is_gzip_path(sketch_path) {
        format!("fasterani-contigs.{content_id}.tsv.bgz")
    } else {
        format!("fasterani-contigs.{content_id}.tsv")
    };
    sketch_path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .map(|parent| parent.join(&filename))
        .unwrap_or_else(|| PathBuf::from(filename))
}

struct Fnv128Writer {
    state: u128,
}

impl Default for Fnv128Writer {
    fn default() -> Self {
        Self {
            state: 0x6c62_272e_07bb_0142_62b8_2175_6295_c58d,
        }
    }
}

impl Write for Fnv128Writer {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        const FNV_PRIME: u128 = 0x0000_0000_0100_0000_0000_0000_0000_013b;
        for byte in bytes {
            self.state ^= u128::from(*byte);
            self.state = self.state.wrapping_mul(FNV_PRIME);
        }
        Ok(bytes.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

fn write_contig_name_records(
    writer: &mut dyn Write,
    files: &[ReferenceFile],
    contig_names: &[ReferenceContigName],
) -> io::Result<()> {
    writeln!(
        writer,
        "contig_id\treference_file_id\treference_file\treference_contig\tsegment_start\tsegment_end"
    )?;

    for (contig_id, contig) in contig_names.iter().enumerate() {
        let reference_file: &str = files
            .get(contig.file_id)
            .map_or("unknown", |file| file.path.as_str());
        writeln!(
            writer,
            "{contig_id}\t{}\t{}\t{}\t{}\t{}",
            contig.file_id,
            tsv_field(reference_file),
            tsv_field(&contig.name),
            contig.segment_start,
            contig.segment_end
        )?;
    }

    Ok(())
}

pub(crate) fn new_contig_sidecar_path(
    sketch_path: &Path,
    files: &[ReferenceFile],
    contig_names: &[ReferenceContigName],
) -> io::Result<PathBuf> {
    let mut digest = Fnv128Writer::default();
    write_contig_name_records(&mut digest, files, contig_names)?;
    Ok(content_addressed_contig_sidecar_path(
        sketch_path,
        &format!("{:032x}", digest.state),
    ))
}

pub(crate) fn contig_sidecar_filename(sidecar_path: &Path) -> io::Result<String> {
    sidecar_path.file_name().map_or_else(
        || {
            Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "contig sidecar path {} has no filename",
                    sidecar_path.display()
                ),
            ))
        },
        |name| Ok(name.to_string_lossy().into_owned()),
    )
}

pub(crate) fn contig_sidecar_entry_path(sketch_path: &Path, filename: &str) -> io::Result<PathBuf> {
    let filename_path: PathBuf = PathBuf::from(filename);
    let mut components = filename_path.components();
    if !matches!(components.next(), Some(Component::Normal(_))) || components.next().is_some() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("invalid contig sidecar filename in sketch metadata: {filename:?}"),
        ));
    }

    Ok(sketch_path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .map(|parent| parent.join(&filename_path))
        .unwrap_or(filename_path))
}

pub(crate) fn manifest_path(prefix: &Path) -> PathBuf {
    append_path_suffix(prefix, ".manifest.json")
}

pub(crate) fn global_frequency_path(prefix: &Path, generation_id: &str) -> PathBuf {
    append_path_suffix(prefix, &format!(".{generation_id}.frequencies.bin"))
}

pub(crate) fn global_frequency_filename(prefix: &Path, generation_id: &str) -> String {
    let path: PathBuf = global_frequency_path(prefix, generation_id);
    path.file_name().map_or_else(
        || path.to_string_lossy().into_owned(),
        |name| name.to_string_lossy().into_owned(),
    )
}

pub(crate) fn global_frequency_entry_path(prefix: &Path, filename: &str) -> PathBuf {
    let filename_path: PathBuf = PathBuf::from(filename);
    if filename_path.is_absolute() {
        return filename_path;
    }

    prefix
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .map(|parent| parent.join(filename))
        .unwrap_or(filename_path)
}

pub(crate) fn shard_path(
    prefix: &Path,
    generation_id: &str,
    shard_index: usize,
    bgzip: bool,
) -> PathBuf {
    if bgzip {
        append_path_suffix(
            prefix,
            &format!(".{generation_id}.{shard_index}.fasketch.bgz"),
        )
    } else {
        append_path_suffix(prefix, &format!(".{generation_id}.{shard_index}.fasketch"))
    }
}

pub(crate) fn shard_filename(
    prefix: &Path,
    generation_id: &str,
    shard_index: usize,
    bgzip: bool,
) -> String {
    shard_path(prefix, generation_id, shard_index, bgzip)
        .file_name()
        .map_or_else(
            || {
                shard_path(prefix, generation_id, shard_index, bgzip)
                    .to_string_lossy()
                    .into_owned()
            },
            |name| name.to_string_lossy().into_owned(),
        )
}

pub(crate) fn shard_entry_path(prefix: &Path, entry: &ShardManifestEntry) -> PathBuf {
    let filename_path: PathBuf = PathBuf::from(&entry.filename);
    if filename_path.is_absolute() {
        return filename_path;
    }

    prefix
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .map(|parent| parent.join(&entry.filename))
        .unwrap_or(filename_path)
}

pub(crate) fn legacy_sketch_path(prefix: &Path) -> Option<PathBuf> {
    if prefix.exists() {
        return Some(prefix.to_path_buf());
    }

    if prefix.extension().is_none() {
        for extension in ["fasketch.bgz", "fasketch"] {
            let path: PathBuf = append_path_suffix(prefix, &format!(".{extension}"));
            if path.exists() {
                return Some(path);
            }
        }
    }

    None
}

pub(crate) fn reference_list_checksum(references: &[FastaInput]) -> u64 {
    const FNV_OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
    const FNV_PRIME: u64 = 0x0100_0000_01b3;

    let mut hash: u64 = FNV_OFFSET;
    for reference in references {
        for byte in reference
            .label
            .as_bytes()
            .iter()
            .copied()
            .chain(std::iter::once(0))
        {
            hash ^= u64::from(byte);
            hash = hash.wrapping_mul(FNV_PRIME);
        }
    }

    hash
}

fn tsv_field(value: &str) -> String {
    value
        .chars()
        .map(|ch| match ch {
            '\t' | '\n' | '\r' => ' ',
            _ => ch,
        })
        .collect()
}

fn files_equal(left: &Path, right: &Path) -> io::Result<bool> {
    if fs::metadata(left)?.len() != fs::metadata(right)?.len() {
        return Ok(false);
    }

    let mut left_reader = BufReader::new(fs::File::open(left)?);
    let mut right_reader = BufReader::new(fs::File::open(right)?);
    let mut left_buffer = [0u8; 64 * 1024];
    let mut right_buffer = [0u8; 64 * 1024];
    loop {
        let left_read: usize = left_reader.read(&mut left_buffer)?;
        let right_read: usize = right_reader.read(&mut right_buffer)?;
        if left_read != right_read || left_buffer[..left_read] != right_buffer[..right_read] {
            return Ok(false);
        }
        if left_read == 0 {
            return Ok(true);
        }
    }
}

fn publish_immutable_temp_file(temp_path: &Path, final_path: &Path) -> io::Result<()> {
    fs::File::open(temp_path)?.sync_all()?;
    match fs::hard_link(temp_path, final_path) {
        Ok(()) => {
            fs::remove_file(temp_path)?;
            sync_parent(final_path)
        }
        Err(err) if err.kind() == io::ErrorKind::AlreadyExists => {
            let identical: bool = files_equal(temp_path, final_path)?;
            fs::remove_file(temp_path)?;
            if identical {
                Ok(())
            } else {
                Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "content-addressed contig sidecar {} already exists with different bytes",
                        final_path.display()
                    ),
                ))
            }
        }
        Err(err) => Err(err),
    }
}

pub(crate) fn write_contig_name_sidecar(
    sidecar_path: &Path,
    files: &[ReferenceFile],
    contig_names: &[ReferenceContigName],
) -> io::Result<u64> {
    let path: PathBuf = sidecar_path.to_path_buf();
    if let Some(parent) = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        fs::create_dir_all(parent)?;
    }

    let (temp_path, file): (PathBuf, fs::File) = create_sibling_temp_file(&path)?;
    let result: io::Result<()> = (|| {
        if is_gzip_path(&path) {
            let mut writer = ZBuilder::<Bgzf, fs::File>::new()
                .num_threads(1)
                .from_writer(file);
            write_contig_name_records(&mut *writer, files, contig_names)?;
            writer.finish().map_err(gzp_error_to_io)?;
        } else {
            let mut writer: BufWriter<fs::File> = BufWriter::new(file);
            write_contig_name_records(&mut writer, files, contig_names)?;
            writer.flush()?;
        }
        publish_immutable_temp_file(&temp_path, &path)
    })();
    if result.is_err() {
        let _ = fs::remove_file(&temp_path);
    }
    result?;
    Ok(fs::metadata(&path)?.len())
}

fn load_contig_sidecar_contents(path: &Path) -> io::Result<String> {
    read_text_maybe_gzip(path).map_err(|err| {
        io::Error::new(
            err.kind(),
            format!(
                "failed to read contig sidecar {}; rebuild the sketch database to create it: {err}",
                path.display()
            ),
        )
    })
}

fn parse_usize_field(value: &str, field: &str, path: &Path) -> io::Result<usize> {
    value.parse::<usize>().map_err(|err| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "failed to parse {field} in contig sidecar {}: {err}",
                path.display()
            ),
        )
    })
}

fn parse_u32_field(value: &str, field: &str, path: &Path) -> io::Result<u32> {
    value.parse::<u32>().map_err(|err| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "failed to parse {field} in contig sidecar {}: {err}",
                path.display()
            ),
        )
    })
}

pub(crate) fn load_contig_name_sidecar(
    sidecar_path: &Path,
    expected_contigs: usize,
    expected_file_bytes: u64,
    expected_files: &[ReferenceFile],
    expected_contig_records: &[ContigRecord],
) -> io::Result<Vec<ReferenceContigName>> {
    let path: PathBuf = sidecar_path.to_path_buf();
    let actual_file_bytes: u64 = fs::metadata(&path)
        .map_err(|err| {
            io::Error::new(
                err.kind(),
                format!(
                    "failed to inspect contig sidecar {}; rebuild the sketch database to create it: {err}",
                    path.display()
                ),
            )
        })?
        .len();
    if actual_file_bytes != expected_file_bytes {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "contig sidecar {} has {actual_file_bytes} bytes, expected {expected_file_bytes}; rebuild the sketch database",
                path.display()
            ),
        ));
    }
    let contents: String = load_contig_sidecar_contents(&path)?;
    let mut digest = Fnv128Writer::default();
    digest.write_all(contents.as_bytes())?;
    let expected_content_path: PathBuf =
        content_addressed_contig_sidecar_path(&path, &format!("{:032x}", digest.state));
    if expected_content_path.file_name() != path.file_name() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "contig sidecar {} does not match its content digest; rebuild the sketch database",
                path.display()
            ),
        ));
    }
    if expected_contig_records.len() != expected_contigs {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "reference sketch has {} binary contig records, expected {expected_contigs}",
                expected_contig_records.len()
            ),
        ));
    }
    let mut lines = contents.lines();
    let header: &str = lines.next().ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("contig sidecar {} is empty", path.display()),
        )
    })?;
    if header
        != "contig_id\treference_file_id\treference_file\treference_contig\tsegment_start\tsegment_end"
    {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("contig sidecar {} has an unexpected header", path.display()),
        ));
    }

    let mut contigs: Vec<ReferenceContigName> = Vec::with_capacity(expected_contigs);
    for (line_index, line) in lines.enumerate() {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 6 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "contig sidecar {} line {} has {} fields, expected 6",
                    path.display(),
                    line_index + 2,
                    fields.len()
                ),
            ));
        }

        let contig_id: usize = parse_usize_field(fields[0], "contig_id", &path)?;
        if contig_id != contigs.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "contig sidecar {} expected contig_id {}, found {contig_id}",
                    path.display(),
                    contigs.len()
                ),
            ));
        }

        let file_id: usize = parse_usize_field(fields[1], "reference_file_id", &path)?;
        let expected_record: &ContigRecord = &expected_contig_records[contig_id];
        if file_id != expected_record.file_id as usize {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "contig sidecar {} contig_id {contig_id} has reference_file_id {file_id}, but the sketch record has {}",
                    path.display(),
                    expected_record.file_id
                ),
            ));
        }
        let expected_file: &ReferenceFile = expected_files.get(file_id).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "contig sidecar {} contig_id {contig_id} references missing file_id {file_id}",
                    path.display()
                ),
            )
        })?;
        let expected_reference_file: String = tsv_field(&expected_file.path);
        if fields[2] != expected_reference_file {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "contig sidecar {} contig_id {contig_id} reference_file does not match sketch file_id {file_id}",
                    path.display()
                ),
            ));
        }
        let segment_start: u32 = parse_u32_field(fields[4], "segment_start", &path)?;
        let segment_end: u32 = parse_u32_field(fields[5], "segment_end", &path)?;
        if segment_start > segment_end {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "contig sidecar {} contig_id {contig_id} has segment_start {segment_start} after segment_end {segment_end}",
                    path.display()
                ),
            ));
        }

        contigs.push(ReferenceContigName {
            file_id,
            name: fields[3].to_string(),
            segment_start,
            segment_end,
        });
    }

    if contigs.len() != expected_contigs {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "contig sidecar {} has {} contigs, expected {expected_contigs}",
                path.display(),
                contigs.len()
            ),
        ));
    }

    Ok(contigs)
}

pub(crate) fn unix_timestamp_seconds() -> io::Result<u64> {
    Ok(SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_err(|err| io::Error::other(format!("system clock is before UNIX epoch: {err}")))?
        .as_secs())
}

pub(crate) fn build_generation_id() -> io::Result<String> {
    let nanos: u128 = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_err(|err| io::Error::other(format!("system clock is before UNIX epoch: {err}")))?
        .as_nanos();
    Ok(format!("{nanos:x}-{}", std::process::id()))
}

#[cfg(test)]
mod tests {
    use super::{write_bytes_atomically, SketchOutput};
    use crate::ani::{manifest_path, shard_filename, shard_path, sketch_reference_name};
    use std::{
        env, fs, io,
        io::Write,
        path::PathBuf,
        time::{SystemTime, UNIX_EPOCH},
    };

    fn unique_test_dir(label: &str) -> PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock after epoch")
            .as_nanos();
        env::temp_dir().join(format!("fasterani-{label}-{}-{nanos}", std::process::id()))
    }

    #[test]
    fn sharded_sketch_paths_use_prefix_suffixes() {
        let prefix: PathBuf = PathBuf::from("/tmp/database");

        assert_eq!(
            manifest_path(&prefix),
            PathBuf::from("/tmp/database.manifest.json")
        );
        assert_eq!(
            shard_path(&prefix, "generation", 2, false),
            PathBuf::from("/tmp/database.generation.2.fasketch")
        );
        assert_eq!(
            shard_filename(&prefix, "generation", 2, false),
            "database.generation.2.fasketch"
        );
        assert_eq!(
            shard_path(&prefix, "generation", 2, true),
            PathBuf::from("/tmp/database.generation.2.fasketch.bgz")
        );
        assert_eq!(
            shard_filename(&prefix, "generation", 2, true),
            "database.generation.2.fasketch.bgz"
        );
    }

    #[test]
    fn sketch_reference_name_keeps_only_basename() {
        assert_eq!(
            sketch_reference_name("/tmp/reference/GCF_000146045.2_R64_genomic.fna"),
            "GCF_000146045.2_R64_genomic.fna"
        );
        assert_eq!(sketch_reference_name("relative.fa"), "relative.fa");
    }

    #[test]
    fn interrupted_atomic_output_preserves_the_published_file() -> io::Result<()> {
        let directory = unique_test_dir("atomic-output-interruption");
        fs::create_dir_all(&directory)?;
        let path = directory.join("artifact.bin");
        fs::write(&path, b"published-generation")?;

        {
            let mut output = SketchOutput::create(&path, None, false, 1)?;
            output.writer_mut()?.write_all(b"incomplete-generation")?;
            // Dropping before finish models a failed build. The sibling temporary
            // file is removed and the previously published artifact is untouched.
        }

        assert_eq!(fs::read(&path)?, b"published-generation");
        assert_eq!(fs::read_dir(&directory)?.count(), 1);
        fs::remove_dir_all(directory)?;
        Ok(())
    }

    #[test]
    fn atomic_byte_write_replaces_only_after_complete_write() -> io::Result<()> {
        let directory = unique_test_dir("atomic-byte-write");
        fs::create_dir_all(&directory)?;
        let path = directory.join("manifest.json");
        fs::write(&path, b"old-manifest")?;

        write_bytes_atomically(&path, b"new-manifest")?;

        assert_eq!(fs::read(&path)?, b"new-manifest");
        assert_eq!(fs::read_dir(&directory)?.count(), 1);
        fs::remove_dir_all(directory)?;
        Ok(())
    }
}
