//! Subcommand policy; the legacy option-only interface remains available.
use super::input::resolve_path_from_base;
use crate::ani::{
    io_util::FastaInput,
    params_file::ParamsFileConfig,
    sketch::serialize::{legacy_sketch_path, manifest_path},
};
use std::{fs, io, path::Path};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum Command {
    Legacy,
    Query,
    Sketch,
    Update,
    Inspect,
}
impl Command {
    pub(crate) fn parse(value: &str) -> io::Result<Self> {
        match value {
            "query" => Ok(Self::Query),
            "sketch" => Ok(Self::Sketch),
            "update" => Ok(Self::Update),
            "inspect" => Ok(Self::Inspect),
            _ => Err(invalid(format!(
                "unknown command {value:?}; use query, sketch, update, or inspect"
            ))),
        }
    }
    pub(crate) fn take(args: &mut Vec<String>) -> io::Result<Option<Self>> {
        for arg in args.iter_mut() {
            if arg == "--params" {
                *arg = "--params-file".into();
            }
        }
        if args.first().is_some_and(|value| {
            matches!(value.as_str(), "query" | "sketch" | "update" | "inspect")
        }) {
            Ok(Some(Self::parse(&args.remove(0))?))
        } else {
            Ok(None)
        }
    }
    pub(crate) fn name(self) -> Option<&'static str> {
        match self {
            Self::Legacy => None,
            Self::Query => Some("query"),
            Self::Sketch => Some("sketch"),
            Self::Update => Some("update"),
            Self::Inspect => Some("inspect"),
        }
    }
}
pub(super) fn invalid(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidInput, message.into())
}
pub(super) fn validate_file_inputs(command: Command, p: &ParamsFileConfig) -> io::Result<()> {
    if matches!(command, Command::Sketch | Command::Update) && p.reference_files.is_some() {
        return Err(invalid(
            "sketch/update require newline-separated list files, not reference_files",
        ));
    }
    if command != Command::Update && (p.add_lists.is_some() || p.remove_lists.is_some()) {
        return Err(invalid("add_lists/remove_lists require the update command"));
    }
    if command == Command::Update && p.reference_lists.is_some() {
        return Err(invalid(
            "update uses add_lists and remove_lists, not reference_lists",
        ));
    }
    Ok(())
}
pub(super) fn validate_flag(command: Command, flag: &str) -> io::Result<()> {
    let bad = match flag {
        "--reference" => matches!(
            command,
            Command::Sketch | Command::Update | Command::Inspect
        ),
        "--reference-list" => matches!(command, Command::Update | Command::Inspect),
        "--add-list" | "--remove-list" => command != Command::Update,
        "--output" => command != Command::Sketch,
        "--query" | "--query-list" | "--query-name" => {
            matches!(command, Command::Sketch | Command::Inspect)
        }
        "--shards" => matches!(
            command,
            Command::Sketch | Command::Update | Command::Inspect
        ),
        _ => false,
    };
    if bad {
        return Err(invalid(format!(
            "{flag} is not supported by {}; use list inputs for sketch/update",
            command.name().unwrap_or("legacy mode")
        )));
    }
    Ok(())
}
pub(super) fn add_removals(
    value: &str,
    base: Option<&Path>,
    ids: &mut Vec<String>,
) -> io::Result<()> {
    let path = resolve_path_from_base(value, base)?;
    let mut seen: std::collections::HashSet<String> = ids.iter().cloned().collect();
    for line in fs::read_to_string(path)?.lines() {
        let id = line.trim();
        if id.is_empty() || id.starts_with('#') {
            continue;
        }
        // Exact stored identifiers, not paths that must still exist.
        if !seen.insert(id.to_owned()) {
            return Err(invalid(format!("duplicate removal identifier: {id}")));
        }
        ids.push(id.to_owned());
    }
    Ok(())
}
pub(super) fn validate_operation(
    command: Command,
    refs: &[FastaInput],
    queries: &[FastaInput],
    prefix: Option<&Path>,
    removals: &[String],
    shards: bool,
) -> io::Result<()> {
    let exists =
        prefix.is_some_and(|p| manifest_path(p).exists() || legacy_sketch_path(p).is_some());
    match command {
        Command::Legacy => {}
        Command::Query => {
            if queries.is_empty() {
                return Err(invalid("query requires --query or --query-list"));
            }
            if prefix.is_some() && !exists {
                return Err(invalid(
                    "saved database does not exist; build it with sketch first",
                ));
            }
        }
        Command::Sketch => {
            if prefix.is_none() || refs.is_empty() {
                return Err(invalid("sketch requires --reference-list and --output"));
            }
            if exists {
                return Err(invalid("database already exists; use update"));
            }
            if !queries.is_empty() {
                return Err(invalid("sketch does not accept queries; use query"));
            }
        }
        Command::Update => {
            if !exists {
                return Err(invalid(
                    "update requires an existing --reference-sketch database",
                ));
            }
            if refs.is_empty() && removals.is_empty() {
                return Err(invalid(
                    "update requires a nonempty --add-list and/or --remove-list",
                ));
            }
        }
        Command::Inspect => {
            if !exists {
                return Err(invalid(
                    "inspect requires an existing --reference-sketch database",
                ));
            }
            if !refs.is_empty() || !queries.is_empty() {
                return Err(invalid("inspect does not accept FASTA inputs"));
            }
        }
    }
    if shards
        && matches!(
            command,
            Command::Sketch | Command::Update | Command::Inspect
        )
    {
        return Err(invalid("--shards is a query option"));
    }
    Ok(())
}
