//! Reading and writing sketch files, shard manifests, and contig sidecars.

use std::{
    fs, io,
    io::{BufWriter, Write},
    path::{Path, PathBuf},
    time::{SystemTime, UNIX_EPOCH},
};

use gzp::{deflate::Bgzf, ZBuilder};

use crate::ani::{
    append_path_suffix, compress_file_to_bgzf, gzp_error_to_io, is_gzip_path, read_text_maybe_gzip,
    ReferenceContigName, ReferenceFile, ScratchFile, ShardManifestEntry,
};

pub(crate) struct SketchOutput {
    pub(crate) final_path: PathBuf,
    pub(crate) write_path: PathBuf,
    pub(crate) scratch: Option<ScratchFile>,
    pub(crate) writer: Option<BufWriter<fs::File>>,
    pub(crate) bgzip: bool,
    pub(crate) threads: usize,
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

        if bgzip {
            let (scratch, file): (ScratchFile, fs::File) =
                ScratchFile::create(tmp_dir, "uncompressed-sketch")?;
            let write_path: PathBuf = scratch.path.clone();
            Ok(Self {
                final_path: final_path.to_path_buf(),
                write_path,
                scratch: Some(scratch),
                writer: Some(BufWriter::new(file)),
                bgzip,
                threads,
            })
        } else {
            let file: fs::File = fs::File::create(final_path)?;
            Ok(Self {
                final_path: final_path.to_path_buf(),
                write_path: final_path.to_path_buf(),
                scratch: None,
                writer: Some(BufWriter::new(file)),
                bgzip,
                threads,
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
            compress_file_to_bgzf(&self.write_path, &self.final_path, self.threads)?;
        }

        let output_len: u64 = fs::metadata(&self.final_path)?.len();
        drop(self.scratch.take());
        Ok(output_len)
    }
}

pub(crate) fn contig_sidecar_path(sketch_path: &Path) -> PathBuf {
    let sketch_path_string: String = sketch_path.to_string_lossy().into_owned();
    if let Some(uncompressed_name) = sketch_path_string.strip_suffix(".bgz") {
        return append_path_suffix(Path::new(uncompressed_name), ".contigs.tsv.bgz");
    }

    append_path_suffix(sketch_path, ".contigs.tsv")
}

pub(crate) fn manifest_path(prefix: &Path) -> PathBuf {
    append_path_suffix(prefix, ".manifest.json")
}

pub(crate) fn shard_path(prefix: &Path, shard_index: usize, bgzip: bool) -> PathBuf {
    if bgzip {
        append_path_suffix(prefix, &format!(".{shard_index}.fasketch.bgz"))
    } else {
        append_path_suffix(prefix, &format!(".{shard_index}.fasketch"))
    }
}

pub(crate) fn shard_filename(prefix: &Path, shard_index: usize, bgzip: bool) -> String {
    shard_path(prefix, shard_index, bgzip)
        .file_name()
        .map(|name| name.to_string_lossy().into_owned())
        .unwrap_or_else(|| {
            shard_path(prefix, shard_index, bgzip)
                .to_string_lossy()
                .into_owned()
        })
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

pub(crate) fn reference_list_checksum(reference_paths: &[String]) -> u64 {
    const FNV_OFFSET: u64 = 0xcbf29ce484222325;
    const FNV_PRIME: u64 = 0x100000001b3;

    let mut hash: u64 = FNV_OFFSET;
    for path in reference_paths {
        for byte in path.as_bytes().iter().copied().chain(std::iter::once(0)) {
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

pub(crate) fn write_contig_name_sidecar(
    sketch_path: &Path,
    files: &[ReferenceFile],
    contig_names: &[ReferenceContigName],
    threads: usize,
) -> io::Result<()> {
    let path: PathBuf = contig_sidecar_path(sketch_path);
    if let Some(parent) = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        fs::create_dir_all(parent)?;
    }

    let write_records = |writer: &mut dyn Write| -> io::Result<()> {
        writeln!(
            writer,
            "contig_id\treference_file_id\treference_file\treference_contig\tsegment_start\tsegment_end"
        )?;

        for (contig_id, contig) in contig_names.iter().enumerate() {
            let reference_file: &str = files
                .get(contig.file_id)
                .map(|file| file.path.as_str())
                .unwrap_or("unknown");
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
    };

    if is_gzip_path(&path) {
        let file: fs::File = fs::File::create(path)?;
        let mut writer = ZBuilder::<Bgzf, fs::File>::new()
            .num_threads(threads.max(1))
            .from_writer(file);
        write_records(&mut *writer)?;
        writer.finish().map_err(gzp_error_to_io)?;
    } else {
        let file: fs::File = fs::File::create(path)?;
        let mut writer: BufWriter<fs::File> = BufWriter::new(file);
        write_records(&mut writer)?;
        writer.flush()?;
    }

    Ok(())
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

fn uncompressed_contig_sidecar_path_for_bgzip_sketch(sketch_path: &Path) -> Option<PathBuf> {
    let sketch_path_string: String = sketch_path.to_string_lossy().into_owned();
    sketch_path_string
        .strip_suffix(".bgz")
        .map(|uncompressed_name| append_path_suffix(Path::new(uncompressed_name), ".contigs.tsv"))
}

fn read_contig_sidecar_text(sketch_path: &Path) -> io::Result<(PathBuf, String)> {
    let path: PathBuf = contig_sidecar_path(sketch_path);
    match load_contig_sidecar_contents(&path) {
        Ok(contents) => Ok((path, contents)),
        Err(error) if error.kind() == io::ErrorKind::NotFound => {
            if let Some(fallback_path) =
                uncompressed_contig_sidecar_path_for_bgzip_sketch(sketch_path)
            {
                let contents: String = load_contig_sidecar_contents(&fallback_path)?;
                Ok((fallback_path, contents))
            } else {
                Err(error)
            }
        }
        Err(error) => Err(error),
    }
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
    sketch_path: &Path,
    expected_contigs: usize,
) -> io::Result<Vec<ReferenceContigName>> {
    let (path, contents): (PathBuf, String) = read_contig_sidecar_text(sketch_path)?;
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

        contigs.push(ReferenceContigName {
            file_id: parse_usize_field(fields[1], "reference_file_id", &path)?,
            name: fields[3].to_string(),
            segment_start: parse_u32_field(fields[4], "segment_start", &path)?,
            segment_end: parse_u32_field(fields[5], "segment_end", &path)?,
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

#[cfg(test)]
mod tests {
    use crate::ani::{manifest_path, shard_filename, shard_path, sketch_reference_name};
    use std::path::PathBuf;

    #[test]
    fn sharded_sketch_paths_use_prefix_suffixes() {
        let prefix: PathBuf = PathBuf::from("/tmp/database");

        assert_eq!(
            manifest_path(&prefix),
            PathBuf::from("/tmp/database.manifest.json")
        );
        assert_eq!(
            shard_path(&prefix, 2, false),
            PathBuf::from("/tmp/database.2.fasketch")
        );
        assert_eq!(shard_filename(&prefix, 2, false), "database.2.fasketch");
        assert_eq!(
            shard_path(&prefix, 2, true),
            PathBuf::from("/tmp/database.2.fasketch.bgz")
        );
        assert_eq!(shard_filename(&prefix, 2, true), "database.2.fasketch.bgz");
    }

    #[test]
    fn sketch_reference_name_keeps_only_basename() {
        assert_eq!(
            sketch_reference_name("/tmp/reference/GCF_000146045.2_R64_genomic.fna"),
            "GCF_000146045.2_R64_genomic.fna"
        );
        assert_eq!(sketch_reference_name("relative.fa"), "relative.fa");
    }
}
