//! FASTA input validation and path/list resolution.

use std::{
    fs, io, path,
    path::{Path, PathBuf},
};

use crate::ani::io_util::{is_stdin_path, FastaInput};

use super::effective_config::{ParameterSource, RuntimeStartupOutput, StartupValue};

const MIN_FASTA_FILE_BYTES: u64 = 100;

fn validate_fasta_file_path(
    path: &Path,
    input_kind: &str,
    original_path: &str,
    list_path: Option<&str>,
) -> io::Result<()> {
    match fs::metadata(path) {
        Ok(meta) if meta.is_file() && meta.len() > MIN_FASTA_FILE_BYTES => Ok(()),
        Ok(meta) if meta.is_file() => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            match list_path {
                Some(_) => format!(
                    "Path in list is too small: {} ({} bytes, must be > {MIN_FASTA_FILE_BYTES})",
                    original_path,
                    meta.len()
                ),
                None => format!(
                    "{input_kind} file is too small: {} ({} bytes, must be > {MIN_FASTA_FILE_BYTES})",
                    path.display(),
                    meta.len()
                ),
            },
        )),
        Ok(_) => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            match list_path {
                Some(_) => format!("Path in list is not a file: {original_path}"),
                None => format!("{input_kind} path is not a file: {}", path.display()),
            },
        )),
        Err(err) => Err(io::Error::new(
            err.kind(),
            match list_path {
                Some(list_path) => format!(
                    "Cannot access path '{original_path}' from list: {list_path}\n{err}"
                ),
                None => format!(
                    "Cannot access {} file {}: {}",
                    input_kind.to_ascii_lowercase(),
                    path.display(),
                    err
                ),
            },
        )),
    }
}

fn validate_and_read_path_list_from_base(
    list_path: &str,
    base_dir: Option<&Path>,
    skip_validation: bool,
) -> io::Result<(PathBuf, Vec<String>)> {
    let absolute_path = resolve_path_from_base(list_path, base_dir)?;
    let contents = fs::read_to_string(&absolute_path)?;
    let mut valid_paths = Vec::new();
    let list_base_dir = if base_dir.is_some() {
        absolute_path.parent()
    } else {
        None
    };

    for line in contents.lines() {
        let path = line.trim();
        if path.is_empty() || path.starts_with('#') {
            continue;
        }

        let absolute_item_path = resolve_path_from_base(path, list_base_dir)?;
        if !skip_validation {
            validate_fasta_file_path(&absolute_item_path, "Input", path, Some(list_path))?;
        }
        valid_paths.push(absolute_item_path.to_string_lossy().into_owned());
    }
    Ok((absolute_path, valid_paths))
}

pub(super) fn resolve_path_from_base(value: &str, base_dir: Option<&Path>) -> io::Result<PathBuf> {
    let path = PathBuf::from(value);
    if path.is_absolute() {
        return path::absolute(path);
    }
    match base_dir {
        Some(base_dir) => path::absolute(base_dir.join(path)),
        None => path::absolute(path),
    }
}

pub(super) fn params_file_base_dir(params_file_path: &Path) -> Option<PathBuf> {
    let absolute_path = path::absolute(params_file_path).ok()?;
    absolute_path.parent().map(Path::to_path_buf)
}

pub(super) fn add_reference_file(
    value: &str,
    source: ParameterSource,
    base_dir: Option<&Path>,
    skip_validation: bool,
    references: &mut Vec<FastaInput>,
    startup_output: &mut RuntimeStartupOutput,
) -> io::Result<()> {
    let reference_absolute_path = resolve_path_from_base(value, base_dir)?;
    if !skip_validation {
        validate_fasta_file_path(&reference_absolute_path, "Reference", value, None)?;
    }

    startup_output.reference_files.push(StartupValue::new(
        reference_absolute_path.to_string_lossy().into_owned(),
        source,
    ));
    let input_path = match source {
        ParameterSource::Cli => value.to_owned(),
        ParameterSource::ParamsFile => reference_absolute_path.display().to_string(),
    };
    references.push(FastaInput::from_path(input_path));
    Ok(())
}

pub(super) fn add_reference_list(
    value: &str,
    source: ParameterSource,
    base_dir: Option<&Path>,
    skip_validation: bool,
    references: &mut Vec<FastaInput>,
    startup_output: &mut RuntimeStartupOutput,
) -> io::Result<()> {
    let (absolute_path, validated_paths) =
        validate_and_read_path_list_from_base(value, base_dir, skip_validation)?;
    match fs::exists(&absolute_path) {
        Ok(true) => {}
        Ok(false) => eprintln!(
            "ERROR\tevent=reference_list_missing\tpath={}",
            absolute_path.display()
        ),
        Err(e) => eprintln!(
            "ERROR\tevent=reference_list_load_failed\tpath={}\terror={}",
            absolute_path.display(),
            e
        ),
    }

    startup_output.reference_lists.push(StartupValue::new(
        absolute_path.to_string_lossy().into_owned(),
        source,
    ));
    references.extend(validated_paths.into_iter().map(FastaInput::from_path));
    Ok(())
}

pub(super) fn add_query_file(
    value: &str,
    source: ParameterSource,
    base_dir: Option<&Path>,
    skip_validation: bool,
    queries: &mut Vec<FastaInput>,
    startup_output: &mut RuntimeStartupOutput,
) -> io::Result<()> {
    if is_stdin_path(value) {
        startup_output
            .query_files
            .push(StartupValue::new(value, source));
        queries.push(FastaInput::from_stdin(None));
        return Ok(());
    }

    let query_absolute_path = resolve_path_from_base(value, base_dir)?;
    if !skip_validation {
        validate_fasta_file_path(&query_absolute_path, "Query", value, None)?;
    }
    startup_output.query_files.push(StartupValue::new(
        query_absolute_path.to_string_lossy().into_owned(),
        source,
    ));
    let input_path = match source {
        ParameterSource::Cli => value.to_owned(),
        ParameterSource::ParamsFile => query_absolute_path.display().to_string(),
    };
    queries.push(FastaInput::from_path(input_path));
    Ok(())
}

pub(super) fn add_query_list(
    value: &str,
    source: ParameterSource,
    base_dir: Option<&Path>,
    skip_validation: bool,
    queries: &mut Vec<FastaInput>,
    startup_output: &mut RuntimeStartupOutput,
) -> io::Result<()> {
    let (absolute_path, validated_paths) =
        validate_and_read_path_list_from_base(value, base_dir, skip_validation)?;
    match fs::exists(&absolute_path) {
        Ok(true) => {}
        Ok(false) => eprintln!(
            "ERROR\tevent=query_list_missing\tpath={}",
            absolute_path.display()
        ),
        Err(e) => eprintln!(
            "ERROR\tevent=query_list_load_failed\tpath={}\terror={}",
            absolute_path.display(),
            e
        ),
    }

    startup_output.query_lists.push(StartupValue::new(
        absolute_path.to_string_lossy().into_owned(),
        source,
    ));
    queries.extend(validated_paths.into_iter().map(FastaInput::from_path));
    Ok(())
}
