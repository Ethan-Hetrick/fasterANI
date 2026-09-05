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

use crate::ani::{
    io_util::{
        append_path_suffix, compress_file_to_bgzf, sketch_reference_name, FastaInput, ScratchFile,
    },
    mmap::MmapFile,
    model::reference::{ReferenceContigName, ReferenceFile, ShardManifestEntry},
};

const NAME_SIDECAR_MAGIC: [u8; 8] = *b"FANINAM\0";
const NAME_SIDECAR_VERSION: u32 = 2;
const NAME_SIDECAR_HEADER_BYTES: usize = 72;
const NAME_SIDECAR_GENOME_RECORD_BYTES: usize = 32;
const NAME_SIDECAR_CONTIG_RECORD_BYTES: usize = 24;

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
fn content_addressed_name_sidecar_path(sketch_path: &Path, content_id: &str) -> PathBuf {
    let filename = format!("fasterani-names.{content_id}.bin");
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

fn push_u32(bytes: &mut Vec<u8>, value: u32) {
    bytes.extend_from_slice(&value.to_le_bytes());
}

fn push_u64(bytes: &mut Vec<u8>, value: u64) {
    bytes.extend_from_slice(&value.to_le_bytes());
}

fn name_blob_record(blob: &mut Vec<u8>, name: &str) -> io::Result<(u64, u32)> {
    let offset = u64::try_from(blob.len()).map_err(|err| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("name sidecar blob offset exceeds u64: {err}"),
        )
    })?;
    let length = u32::try_from(name.len()).map_err(|err| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("name sidecar entry exceeds u32 bytes: {err}"),
        )
    })?;
    blob.extend_from_slice(name.as_bytes());
    Ok((offset, length))
}

fn build_name_sidecar_bytes(
    files: &[ReferenceFile],
    contig_names: &[ReferenceContigName],
) -> io::Result<Vec<u8>> {
    let mut blob = Vec::new();
    let mut genome_records = Vec::with_capacity(files.len());
    for file in files {
        let saved_name = sketch_reference_name(&file.path);
        let (name_offset, name_length) = name_blob_record(&mut blob, &saved_name)?;
        genome_records.push((
            name_offset,
            name_length,
            file.original_length,
            file.mapped_length,
        ));
    }
    let mut contig_records = Vec::with_capacity(contig_names.len());
    for contig in contig_names {
        let (name_offset, name_length) = name_blob_record(&mut blob, &contig.name)?;
        let file_id = u32::try_from(contig.file_id).map_err(|err| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("contig file id exceeds u32: {err}"),
            )
        })?;
        contig_records.push((
            name_offset,
            name_length,
            file_id,
            contig.segment_start,
            contig.segment_end,
        ));
    }

    let genome_index_offset = NAME_SIDECAR_HEADER_BYTES;
    let contig_index_offset = genome_index_offset
        .checked_add(
            genome_records
                .len()
                .checked_mul(NAME_SIDECAR_GENOME_RECORD_BYTES)
                .ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        "name sidecar genome index overflow",
                    )
                })?,
        )
        .ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "name sidecar offset overflow")
        })?;
    let blob_offset = contig_index_offset
        .checked_add(
            contig_records
                .len()
                .checked_mul(NAME_SIDECAR_CONTIG_RECORD_BYTES)
                .ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        "name sidecar contig index overflow",
                    )
                })?,
        )
        .ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "name sidecar offset overflow")
        })?;
    let file_size = blob_offset
        .checked_add(blob.len())
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "name sidecar size overflow"))?;

    let mut bytes = Vec::with_capacity(file_size);
    bytes.extend_from_slice(&NAME_SIDECAR_MAGIC);
    push_u32(&mut bytes, NAME_SIDECAR_VERSION);
    push_u32(&mut bytes, NAME_SIDECAR_HEADER_BYTES as u32);
    push_u64(&mut bytes, genome_records.len() as u64);
    push_u64(&mut bytes, contig_records.len() as u64);
    push_u64(&mut bytes, genome_index_offset as u64);
    push_u64(&mut bytes, contig_index_offset as u64);
    push_u64(&mut bytes, blob_offset as u64);
    push_u64(&mut bytes, blob.len() as u64);
    push_u64(&mut bytes, file_size as u64);
    debug_assert_eq!(bytes.len(), NAME_SIDECAR_HEADER_BYTES);
    for (name_offset, name_length, original_length, mapped_length) in genome_records {
        push_u64(&mut bytes, name_offset);
        push_u32(&mut bytes, name_length);
        push_u32(&mut bytes, 0);
        push_u64(&mut bytes, original_length);
        push_u64(&mut bytes, mapped_length);
    }
    for (name_offset, name_length, file_id, segment_start, segment_end) in contig_records {
        push_u64(&mut bytes, name_offset);
        push_u32(&mut bytes, name_length);
        push_u32(&mut bytes, file_id);
        push_u32(&mut bytes, segment_start);
        push_u32(&mut bytes, segment_end);
    }
    bytes.extend_from_slice(&blob);
    Ok(bytes)
}

pub(crate) fn write_name_sidecar(
    sketch_path: &Path,
    files: &[ReferenceFile],
    contig_names: &[ReferenceContigName],
) -> io::Result<(String, u64)> {
    let bytes = build_name_sidecar_bytes(files, contig_names)?;
    let mut digest = Fnv128Writer::default();
    digest.write_all(&bytes)?;
    let path = content_addressed_name_sidecar_path(sketch_path, &format!("{:032x}", digest.state));
    let (temp_path, mut file) = create_sibling_temp_file(&path)?;
    let result = (|| {
        file.write_all(&bytes)?;
        file.flush()?;
        drop(file);
        publish_immutable_temp_file(&temp_path, &path)
    })();
    if result.is_err() {
        let _ = fs::remove_file(&temp_path);
    }
    result?;
    let filename = path
        .file_name()
        .map(|name| name.to_string_lossy().into_owned())
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                "name sidecar path has no filename",
            )
        })?;
    Ok((filename, bytes.len() as u64))
}

#[cfg_attr(not(test), allow(dead_code))]
fn read_u32(bytes: &[u8], offset: usize) -> io::Result<u32> {
    let raw: [u8; 4] = bytes
        .get(offset..offset + 4)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "truncated name sidecar"))?
        .try_into()
        .expect("four-byte slice");
    Ok(u32::from_le_bytes(raw))
}

#[cfg_attr(not(test), allow(dead_code))]
fn read_u64(bytes: &[u8], offset: usize) -> io::Result<u64> {
    let raw: [u8; 8] = bytes
        .get(offset..offset + 8)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "truncated name sidecar"))?
        .try_into()
        .expect("eight-byte slice");
    Ok(u64::from_le_bytes(raw))
}

/// Experimental allocation-free reader for an indexed reference-name sidecar.
#[cfg_attr(not(test), allow(dead_code))]
pub(crate) struct NameSidecar {
    mmap: MmapFile,
    genome_count: usize,
    contig_count: usize,
    genome_index_offset: usize,
    contig_index_offset: usize,
    blob_offset: usize,
    blob_len: usize,
}

#[cfg_attr(not(test), allow(dead_code))]
impl NameSidecar {
    pub(crate) fn open(path: &Path, expected_file_bytes: u64) -> io::Result<Self> {
        let mmap = MmapFile::open(path)?;
        let bytes = mmap.as_slice();
        if expected_file_bytes != 0 && expected_file_bytes != bytes.len() as u64 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "name sidecar byte count mismatch",
            ));
        }
        if bytes.get(..8) != Some(NAME_SIDECAR_MAGIC.as_slice()) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid name sidecar magic",
            ));
        }
        if read_u32(bytes, 8)? != NAME_SIDECAR_VERSION
            || read_u32(bytes, 12)? as usize != NAME_SIDECAR_HEADER_BYTES
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "unsupported name sidecar format",
            ));
        }
        let genome_count = usize::try_from(read_u64(bytes, 16)?).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "name sidecar genome count exceeds usize",
            )
        })?;
        let contig_count = usize::try_from(read_u64(bytes, 24)?).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "name sidecar contig count exceeds usize",
            )
        })?;
        let genome_index_offset = usize::try_from(read_u64(bytes, 32)?).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "name sidecar offset exceeds usize",
            )
        })?;
        let contig_index_offset = usize::try_from(read_u64(bytes, 40)?).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "name sidecar offset exceeds usize",
            )
        })?;
        let blob_offset = usize::try_from(read_u64(bytes, 48)?).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "name sidecar offset exceeds usize",
            )
        })?;
        let blob_len = usize::try_from(read_u64(bytes, 56)?).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "name sidecar blob length exceeds usize",
            )
        })?;
        let declared_size = usize::try_from(read_u64(bytes, 64)?).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "name sidecar size exceeds usize",
            )
        })?;
        let expected_contig_offset = NAME_SIDECAR_HEADER_BYTES
            .checked_add(
                genome_count
                    .checked_mul(NAME_SIDECAR_GENOME_RECORD_BYTES)
                    .ok_or_else(|| {
                        io::Error::new(io::ErrorKind::InvalidData, "name sidecar index overflow")
                    })?,
            )
            .ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "name sidecar index overflow")
            })?;
        let expected_blob_offset = expected_contig_offset
            .checked_add(
                contig_count
                    .checked_mul(NAME_SIDECAR_CONTIG_RECORD_BYTES)
                    .ok_or_else(|| {
                        io::Error::new(io::ErrorKind::InvalidData, "name sidecar index overflow")
                    })?,
            )
            .ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "name sidecar index overflow")
            })?;
        if genome_index_offset != NAME_SIDECAR_HEADER_BYTES
            || contig_index_offset != expected_contig_offset
            || blob_offset != expected_blob_offset
            || blob_offset.checked_add(blob_len) != Some(bytes.len())
            || declared_size != bytes.len()
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid name sidecar section layout",
            ));
        }
        let reader = Self {
            mmap,
            genome_count,
            contig_count,
            genome_index_offset,
            contig_index_offset,
            blob_offset,
            blob_len,
        };
        for id in 0..reader.genome_count {
            reader.genome_name(id)?;
        }
        for id in 0..reader.contig_count {
            reader.contig_name(id)?;
            if reader.contig_file_id(id)? >= reader.genome_count {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "name sidecar contig references a missing genome",
                ));
            }
            let (start, end) = reader.contig_segment(id)?;
            if start > end {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "name sidecar contig segment is inverted",
                ));
            }
        }
        let mut digest = Fnv128Writer::default();
        digest.write_all(reader.mmap.as_slice())?;
        let expected_path =
            content_addressed_name_sidecar_path(path, &format!("{:032x}", digest.state));
        if expected_path.file_name() != path.file_name() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "name sidecar does not match its content digest",
            ));
        }
        Ok(reader)
    }

    pub(crate) fn genome_count(&self) -> usize {
        self.genome_count
    }
    pub(crate) fn contig_count(&self) -> usize {
        self.contig_count
    }

    fn name_at(&self, record_offset: usize) -> io::Result<&str> {
        let offset = usize::try_from(read_u64(self.mmap.as_slice(), record_offset)?)
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "name offset exceeds usize"))?;
        let len = read_u32(self.mmap.as_slice(), record_offset + 8)? as usize;
        let end = offset
            .checked_add(len)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "name range overflow"))?;
        if end > self.blob_len {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "name extends past sidecar blob",
            ));
        }
        std::str::from_utf8(
            &self.mmap.as_slice()[self.blob_offset + offset..self.blob_offset + end],
        )
        .map_err(|err| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid UTF-8 in name sidecar: {err}"),
            )
        })
    }

    pub(crate) fn genome_name(&self, file_id: usize) -> io::Result<&str> {
        if file_id >= self.genome_count {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "genome id out of range",
            ));
        }
        self.name_at(self.genome_index_offset + file_id * NAME_SIDECAR_GENOME_RECORD_BYTES)
    }

    pub(crate) fn genome_length(&self, file_id: usize) -> io::Result<u64> {
        if file_id >= self.genome_count {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "genome id out of range",
            ));
        }
        read_u64(
            self.mmap.as_slice(),
            self.genome_index_offset + file_id * NAME_SIDECAR_GENOME_RECORD_BYTES + 16,
        )
    }

    pub(crate) fn genome_mapped_length(&self, file_id: usize) -> io::Result<u64> {
        if file_id >= self.genome_count {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "genome id out of range",
            ));
        }
        read_u64(
            self.mmap.as_slice(),
            self.genome_index_offset + file_id * NAME_SIDECAR_GENOME_RECORD_BYTES + 24,
        )
    }

    pub(crate) fn contig_name(&self, contig_id: usize) -> io::Result<&str> {
        if contig_id >= self.contig_count {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "contig id out of range",
            ));
        }
        self.name_at(self.contig_index_offset + contig_id * NAME_SIDECAR_CONTIG_RECORD_BYTES)
    }

    pub(crate) fn contig_file_id(&self, contig_id: usize) -> io::Result<usize> {
        if contig_id >= self.contig_count {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "contig id out of range",
            ));
        }
        Ok(read_u32(
            self.mmap.as_slice(),
            self.contig_index_offset + contig_id * NAME_SIDECAR_CONTIG_RECORD_BYTES + 12,
        )? as usize)
    }

    pub(crate) fn contig_segment(&self, contig_id: usize) -> io::Result<(u32, u32)> {
        if contig_id >= self.contig_count {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "contig id out of range",
            ));
        }
        let offset = self.contig_index_offset + contig_id * NAME_SIDECAR_CONTIG_RECORD_BYTES;
        Ok((
            read_u32(self.mmap.as_slice(), offset + 16)?,
            read_u32(self.mmap.as_slice(), offset + 20)?,
        ))
    }
}

pub(crate) fn sidecar_entry_path(sketch_path: &Path, filename: &str) -> io::Result<PathBuf> {
    let filename_path = PathBuf::from(filename);
    let mut components = filename_path.components();
    if !matches!(components.next(), Some(Component::Normal(_))) || components.next().is_some() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("invalid sidecar filename in sketch metadata: {filename:?}"),
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
            .output_label
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
                        "content-addressed sidecar {} already exists with different bytes",
                        final_path.display()
                    ),
                ))
            }
        }
        Err(err) => Err(err),
    }
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
    use super::{
        build_name_sidecar_bytes, write_bytes_atomically, write_name_sidecar, NameSidecar,
        SketchOutput,
    };
    use crate::ani::{
        io_util::sketch_reference_name,
        model::reference::{ReferenceContigName, ReferenceFile},
        sketch::serialize::{manifest_path, shard_filename, shard_path},
    };
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

    fn name_fixture() -> (Vec<ReferenceFile>, Vec<ReferenceContigName>) {
        (
            vec![
                ReferenceFile {
                    path: "alpha.fa".to_string(),
                    mapped_length: 80,
                    original_length: 101,
                },
                ReferenceFile {
                    path: "βeta.fa".to_string(),
                    mapped_length: 190,
                    original_length: 203,
                },
            ],
            vec![
                ReferenceContigName {
                    file_id: 0,
                    name: "".to_string(),
                    segment_start: 0,
                    segment_end: 10,
                },
                ReferenceContigName {
                    file_id: 0,
                    name: "split-name".to_string(),
                    segment_start: 20,
                    segment_end: 30,
                },
                ReferenceContigName {
                    file_id: 1,
                    name: "β-contig".to_string(),
                    segment_start: 0,
                    segment_end: 40,
                },
            ],
        )
    }

    #[test]
    fn indexed_name_sidecar_round_trips_without_owned_names() -> io::Result<()> {
        let directory = unique_test_dir("name-sidecar-round-trip");
        fs::create_dir_all(&directory)?;
        let sketch_path = directory.join("reference.fasketch.gz");
        let (files, contigs) = name_fixture();
        let (filename, file_bytes) = write_name_sidecar(&sketch_path, &files, &contigs)?;
        let path = directory.join(filename);
        assert_eq!(
            path.extension().and_then(|value| value.to_str()),
            Some("bin")
        );

        let sidecar = NameSidecar::open(&path, file_bytes)?;
        assert!(NameSidecar::open(&path, file_bytes + 1).is_err());
        assert_eq!(sidecar.genome_count(), 2);
        assert_eq!(sidecar.contig_count(), 3);
        assert_eq!(sidecar.genome_name(0)?, "alpha.fa");
        assert_eq!(sidecar.genome_name(1)?, "βeta.fa");
        assert_eq!(sidecar.genome_length(1)?, 203);
        assert_eq!(sidecar.genome_mapped_length(1)?, 190);
        assert_eq!(sidecar.contig_name(0)?, "");
        assert_eq!(sidecar.contig_name(2)?, "β-contig");
        let name_ptr = sidecar.genome_name(0)?.as_ptr() as usize;
        let mmap_start = sidecar.mmap.as_slice().as_ptr() as usize;
        assert!((mmap_start..mmap_start + sidecar.mmap.as_slice().len()).contains(&name_ptr));
        assert!(sidecar.genome_name(2).is_err());
        assert!(sidecar.contig_name(3).is_err());
        drop(sidecar);
        fs::remove_dir_all(directory)?;
        Ok(())
    }

    #[test]
    fn indexed_name_sidecar_rejects_corruption() -> io::Result<()> {
        let directory = unique_test_dir("name-sidecar-corruption");
        fs::create_dir_all(&directory)?;
        let (files, contigs) = name_fixture();
        let valid = build_name_sidecar_bytes(&files, &contigs)?;
        for (label, mut bytes) in [
            ("magic", {
                let mut value = valid.clone();
                value[0] ^= 1;
                value
            }),
            ("offset", {
                let mut value = valid.clone();
                value[48..56].copy_from_slice(&0u64.to_le_bytes());
                value
            }),
            ("utf8", {
                let mut value = valid.clone();
                let last = value.len() - 1;
                value[last] = 0xff;
                value
            }),
            ("truncated", valid[..valid.len() - 1].to_vec()),
        ] {
            let path = directory.join(format!("fasterani-names.{label}.bin"));
            fs::write(&path, &bytes)?;
            assert!(
                NameSidecar::open(&path, bytes.len() as u64).is_err(),
                "{label}"
            );
            bytes.clear();
        }
        fs::remove_dir_all(directory)?;
        Ok(())
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
