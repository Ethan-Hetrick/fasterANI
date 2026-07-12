//! Filesystem, gzip/bgzf, byte-slice, and scratch-file helpers.

use std::{
    env, fs, io,
    io::{BufReader, BufWriter, Cursor, Read, Write},
    path::{Path, PathBuf},
    time::{SystemTime, UNIX_EPOCH},
};

use flate2::read::MultiGzDecoder;
use gzp::{deflate::Bgzf, ZBuilder};
use noodles::fasta;

pub(crate) fn align_up(value: usize, alignment: usize) -> usize {
    debug_assert!(alignment.is_power_of_two());
    (value + alignment - 1) & !(alignment - 1)
}

pub(crate) fn checked_section_end(
    offset: usize,
    count: usize,
    item_size: usize,
) -> io::Result<usize> {
    let byte_len = count.checked_mul(item_size).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            "reference sketch section size overflow",
        )
    })?;
    offset.checked_add(byte_len).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            "reference sketch section offset overflow",
        )
    })
}

#[cfg(test)]
pub(crate) fn sketch_reference_name(path: &str) -> String {
    Path::new(path)
        .file_name()
        .map(|name| name.to_string_lossy().into_owned())
        .unwrap_or_else(|| path.to_string())
}

pub(crate) fn append_path_suffix(path: &Path, suffix: &str) -> PathBuf {
    let mut value = path.as_os_str().to_os_string();
    value.push(suffix);
    PathBuf::from(value)
}

fn path_has_extension(path: &Path, expected: &str) -> bool {
    path.extension()
        .and_then(|extension| extension.to_str())
        .is_some_and(|extension| extension.eq_ignore_ascii_case(expected))
}

pub(crate) fn is_gzip_path(path: &Path) -> bool {
    path_has_extension(path, "gz") || path_has_extension(path, "bgz")
}

pub(crate) fn is_stdin_path(path: &str) -> bool {
    path == "-" || path == "/dev/stdin"
}

/// A FASTA input: an on-disk path, or `-` / `/dev/stdin` for a streamed reader.
#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct FastaInput {
    pub(crate) open: String,
    pub(crate) label: String,
}

impl FastaInput {
    pub(crate) fn from_path(path: String) -> Self {
        Self {
            label: path.clone(),
            open: path,
        }
    }

    pub(crate) fn from_stdin(label: Option<String>) -> Self {
        Self {
            open: "-".to_string(),
            label: label.unwrap_or_else(|| "-".to_string()),
        }
    }
}

fn buf_reader_maybe_gzip(mut reader: impl Read + 'static) -> io::Result<Box<dyn io::BufRead>> {
    let mut header: [u8; 2] = [0; 2];
    let header_len: usize = reader.read(&mut header)?;
    let chained: Box<dyn Read> = Box::new(Cursor::new(header[..header_len].to_vec()).chain(reader));
    if header_len >= 2 && header[0] == 0x1f && header[1] == 0x8b {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(chained))))
    } else {
        Ok(Box::new(BufReader::new(chained)))
    }
}

pub(crate) fn gzp_error_to_io(error: gzp::GzpError) -> io::Error {
    io::Error::other(format!("failed to finish BGZF compression: {error}"))
}

pub(crate) fn open_fasta_reader(path: &str) -> io::Result<fasta::io::Reader<Box<dyn io::BufRead>>> {
    let reader: Box<dyn io::BufRead> = if is_stdin_path(path) {
        buf_reader_maybe_gzip(io::stdin().lock())?
    } else {
        let path_ref: &Path = Path::new(path);
        let file: fs::File = fs::File::open(path_ref)?;
        if is_gzip_path(path_ref) {
            Box::new(BufReader::new(MultiGzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        }
    };

    fasta::io::reader::Builder.build_from_reader(reader)
}

pub(crate) fn read_text_maybe_gzip(path: &Path) -> io::Result<String> {
    let file: fs::File = fs::File::open(path)?;
    let mut contents: String = String::new();

    if is_gzip_path(path) {
        let mut reader: BufReader<MultiGzDecoder<fs::File>> =
            BufReader::new(MultiGzDecoder::new(file));
        reader.read_to_string(&mut contents)?;
    } else {
        let mut reader: BufReader<fs::File> = BufReader::new(file);
        reader.read_to_string(&mut contents)?;
    }

    Ok(contents)
}

pub(crate) fn compress_file_to_bgzf(
    source: &Path,
    destination: &Path,
    threads: usize,
) -> io::Result<()> {
    if let Some(parent) = destination
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        fs::create_dir_all(parent)?;
    }

    let input: fs::File = fs::File::open(source)?;
    let output: fs::File = fs::File::create(destination)?;
    let mut reader: BufReader<fs::File> = BufReader::new(input);
    let mut writer = ZBuilder::<Bgzf, fs::File>::new()
        .num_threads(threads.max(1))
        .from_writer(output);
    io::copy(&mut reader, &mut writer)?;
    writer.finish().map_err(gzp_error_to_io)?;

    Ok(())
}

pub(crate) fn decompress_to_scratch(
    source: &Path,
    tmp_dir: Option<&Path>,
    purpose: &str,
) -> io::Result<ScratchFile> {
    let (scratch, scratch_file): (ScratchFile, fs::File) = ScratchFile::create(tmp_dir, purpose)?;
    let input: fs::File = fs::File::open(source)?;
    let mut reader: BufReader<MultiGzDecoder<fs::File>> =
        BufReader::new(MultiGzDecoder::new(input));
    let mut writer: BufWriter<fs::File> = BufWriter::new(scratch_file);
    io::copy(&mut reader, &mut writer)?;
    writer.flush()?;
    drop(writer);

    Ok(scratch)
}

pub(crate) fn slice_as_bytes<T>(slice: &[T]) -> &[u8] {
    unsafe { std::slice::from_raw_parts(slice.as_ptr().cast::<u8>(), std::mem::size_of_val(slice)) }
}

pub(crate) fn slice_as_bytes_mut<T>(slice: &mut [T]) -> &mut [u8] {
    unsafe {
        std::slice::from_raw_parts_mut(
            slice.as_mut_ptr().cast::<u8>(),
            std::mem::size_of_val(slice),
        )
    }
}

pub(crate) fn write_padding(writer: &mut impl Write, len: usize) -> io::Result<()> {
    const ZEROES: [u8; 8] = [0; 8];
    writer.write_all(&ZEROES[..len])
}

pub(crate) struct ScratchFile {
    pub(crate) path: PathBuf,
}

impl ScratchFile {
    pub(crate) fn create(tmp_dir: Option<&Path>, purpose: &str) -> io::Result<(Self, fs::File)> {
        let directory: PathBuf = tmp_dir.map_or_else(env::temp_dir, Path::to_path_buf);
        fs::create_dir_all(&directory)?;

        for attempt in 0..100u32 {
            let timestamp: u128 = SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .map_err(|err| {
                    io::Error::other(format!("system clock is before UNIX epoch: {err}"))
                })?
                .as_nanos();
            let path: PathBuf = directory.join(format!(
                "fasterani-{}-{timestamp}-{purpose}-{attempt}.tmp",
                std::process::id()
            ));
            match fs::OpenOptions::new()
                .create_new(true)
                .read(true)
                .write(true)
                .open(&path)
            {
                Ok(file) => return Ok((Self { path }, file)),
                Err(error) if error.kind() == io::ErrorKind::AlreadyExists => continue,
                Err(error) => return Err(error),
            }
        }

        Err(io::Error::new(
            io::ErrorKind::AlreadyExists,
            "failed to create a unique fasterANI scratch file",
        ))
    }
}

impl Drop for ScratchFile {
    fn drop(&mut self) {
        let _ = fs::remove_file(&self.path);
    }
}
