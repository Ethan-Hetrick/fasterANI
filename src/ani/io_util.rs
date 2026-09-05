//! Filesystem, gzip, byte-slice, and scratch-file helpers.

use std::{
    env, fs, io,
    io::{BufReader, Cursor, Read, Write},
    path::{Path, PathBuf},
    time::{SystemTime, UNIX_EPOCH},
};

use flate2::read::MultiGzDecoder;
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
    pub(crate) input_path: String,
    pub(crate) output_label: String,
}

impl FastaInput {
    pub(crate) fn from_path(path: String) -> Self {
        Self {
            output_label: path.clone(),
            input_path: path,
        }
    }

    pub(crate) fn from_stdin(label: Option<String>) -> Self {
        Self {
            input_path: "-".to_string(),
            output_label: label.unwrap_or_else(|| "-".to_string()),
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
                Err(error) if error.kind() == io::ErrorKind::AlreadyExists => {}
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
