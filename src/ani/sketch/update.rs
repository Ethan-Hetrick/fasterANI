//! Transactional database updates. Existing generation artifacts are immutable.
use super::{
    database::SketchDatabase,
    frequency::{update_global_frequency_artifact, GlobalFrequencyIndex},
    serialize::{
        build_generation_id, manifest_path, reference_list_checksum, shard_entry_path,
        shard_filename, shard_path, unix_timestamp_seconds, write_bytes_atomically,
    },
};
use crate::ani::{
    io_util::{append_path_suffix, sketch_reference_name, FastaInput},
    model::reference::{ReferenceSketch, ShardManifest, ShardedBuildOptions, SketchParams},
    runtime::RuntimeOptions,
};
use std::{collections::HashSet, env, fs, io, io::Write, path::Path};

/// Keep the inode in place: unlinking a lock file allows competing writers to
/// acquire different inodes. The OS releases this advisory lock on process exit.
pub(crate) struct DatabaseWriteLock {
    _file: fs::File,
}
impl DatabaseWriteLock {
    pub(crate) fn acquire(prefix: &Path) -> io::Result<Self> {
        if let Some(parent) = prefix.parent().filter(|p| !p.as_os_str().is_empty()) {
            fs::create_dir_all(parent)?;
        }
        let file = fs::OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .truncate(false)
            .open(append_path_suffix(prefix, ".lock"))?;
        file.try_lock()
            .map_err(|e| io::Error::other(format!("database writer lock unavailable: {e}")))?;
        Ok(Self { _file: file })
    }
}
fn invalid(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidInput, message.into())
}
pub(crate) fn addition_ids(references: &[FastaInput]) -> io::Result<Vec<String>> {
    let mut seen = HashSet::new();
    references
        .iter()
        .map(|reference| {
            let id = sketch_reference_name(&reference.output_label);
            if !seen.insert(id.clone()) {
                return Err(invalid(format!(
                    "duplicate reference identifier {id:?}; identifiers are stored FASTA basenames"
                )));
            }
            Ok(id)
        })
        .collect()
}

impl SketchDatabase {
    pub(crate) fn update(
        prefix: &Path,
        additions: &[FastaInput],
        removals: &[String],
        params: SketchParams,
        options: ShardedBuildOptions<'_>,
        runtime: RuntimeOptions,
    ) -> io::Result<ShardManifest> {
        let _lock = DatabaseWriteLock::acquire(prefix)?;
        let old_bytes = fs::read(manifest_path(prefix))?;
        let old = Self::load_manifest(prefix, params)?;
        let old_frequencies = GlobalFrequencyIndex::load(prefix, &old)?;
        let add_ids = addition_ids(additions)?;
        let remove_set: HashSet<_> = removals.iter().cloned().collect();
        if remove_set.len() != removals.len() {
            return Err(invalid("duplicate removal identifiers"));
        }
        if additions.is_empty() && removals.is_empty() {
            return Err(invalid("empty update"));
        }
        // Validate the complete operation before writing replacement artifacts.
        let mut existing = HashSet::new();
        let mut masks = Vec::new();
        let mut final_ids = Vec::new();
        for entry in &old.shards {
            let sketch =
                ReferenceSketch::load(&shard_entry_path(prefix, entry), params, false, runtime)?;
            let mut keep = Vec::with_capacity(sketch.files.len());
            for file in &sketch.files {
                if !existing.insert(file.path.clone()) {
                    return Err(invalid(format!("ambiguous duplicate stored identifier {:?}; update requires unique reference identifiers", file.path)));
                }
                let retain = !remove_set.contains(&file.path);
                keep.push(retain);
                if retain {
                    final_ids.push(file.path.clone());
                }
            }
            masks.push(keep);
        }
        for id in removals {
            if !existing.contains(id) {
                return Err(invalid(format!(
                    "unknown removal identifier {id:?}; use inspect"
                )));
            }
        }
        for id in &add_ids {
            if existing.contains(id) && !remove_set.contains(id) {
                return Err(invalid(format!("reference {id:?} already exists; explicitly remove it in the same update to replace it")));
            }
        }
        final_ids.extend(add_ids);
        let generation = build_generation_id()?;
        let mut shards = Vec::new();
        let mut removed_shards = Vec::new();
        let mut added_shards = Vec::new();
        for (entry, keep) in old.shards.iter().zip(masks) {
            if keep.iter().all(|&retain| retain) {
                shards.push(entry.clone());
                continue;
            }
            removed_shards.push(entry.clone());
            if !keep.iter().any(|&retain| retain) {
                continue;
            }
            let sketch =
                ReferenceSketch::load(&shard_entry_path(prefix, entry), params, true, runtime)?;
            let path = shard_path(prefix, &generation, entry.shard_index);
            let stats =
                sketch.save_retaining_files(&keep, params, &path, options.tmp_dir, runtime)?;
            let mut replacement = entry.clone();
            replacement.filename = shard_filename(prefix, &generation, entry.shard_index);
            replacement.reference_count = stats.reference_count;
            replacement.reference_contigs = stats.reference_contig_count;
            replacement.mapped_reference_length = stats.mapped_reference_length;
            replacement.reference_minimizers = stats.reference_minimizer_count;
            replacement.unique_minimizers = stats.unique_minimizer_count;
            replacement.file_bytes = fs::metadata(path)?.len();
            replacement.estimated_file_bytes = replacement.file_bytes;
            added_shards.push(replacement.clone());
            shards.push(replacement);
        }
        if !additions.is_empty() {
            let first = old
                .shards
                .iter()
                .map(|s| s.shard_index)
                .max()
                .unwrap_or(0)
                .checked_add(1)
                .ok_or_else(|| invalid("shard index overflow"))?;
            let new = Self::build_shard_files(
                additions,
                params,
                prefix,
                options,
                runtime,
                &generation,
                first,
            )?;
            added_shards.extend(new.iter().cloned());
            shards.extend(new);
        }
        let frequencies = update_global_frequency_artifact(
            prefix,
            &generation,
            &old_frequencies,
            &removed_shards,
            &added_shards,
            params,
            options.tmp_dir,
            runtime,
        )?;
        let mut next = old.clone();
        let mut offset = 0usize;
        for shard in &mut shards {
            shard.first_reference = offset;
            offset = offset
                .checked_add(shard.reference_count)
                .ok_or_else(|| invalid("reference count overflow"))?;
        }
        if offset != final_ids.len() {
            return Err(io::Error::other("updated reference count mismatch"));
        }
        next.total_references = offset;
        next.total_reference_contigs = shards.iter().map(|s| s.reference_contigs).sum();
        next.total_reference_minimizers = shards.iter().map(|s| s.reference_minimizers).sum();
        next.total_mapped_reference_length = shards.iter().map(|s| s.mapped_reference_length).sum();
        next.total_shard_unique_minimizers = shards.iter().map(|s| s.unique_minimizers).sum();
        next.reference_list_checksum = reference_list_checksum(
            &final_ids
                .iter()
                .cloned()
                .map(FastaInput::from_path)
                .collect::<Vec<_>>(),
        );
        next.reference_identifiers = Some(final_ids);
        next.total_unique_minimizers = frequencies.unique_minimizers;
        next.global_frequency_filename = frequencies.filename;
        next.global_frequency_file_bytes = frequencies.file_bytes;
        next.generation_id = generation;
        next.build_args = env::args().collect();
        next.build_unix_seconds = unix_timestamp_seconds()?;
        next.max_shard_size_bytes = options.max_shard_size_bytes;
        next.shards = shards;
        // Validate the new frequency artifact before publishing the manifest.
        GlobalFrequencyIndex::load(prefix, &next)?;
        // The lock serializes our writers. Also catch external manifest changes.
        if fs::read(manifest_path(prefix))? != old_bytes {
            return Err(io::Error::other("database changed during update"));
        }
        let snapshot = append_path_suffix(prefix, &format!(".{}.manifest.json", old.generation_id));
        if snapshot.exists() {
            if fs::read(&snapshot)? != old_bytes {
                return Err(io::Error::other(
                    "existing generation snapshot differs from current manifest",
                ));
            }
        } else {
            write_bytes_atomically(&snapshot, &old_bytes)?;
        }
        Self::write_manifest(prefix, &next)?;
        Ok(next)
    }

    pub(crate) fn inspect(
        prefix: &Path,
        params: SketchParams,
        runtime: RuntimeOptions,
        writer: &mut dyn Write,
    ) -> io::Result<()> {
        let manifest = Self::load_manifest(prefix, params)?;
        writeln!(writer, "# generation={} k={} w={} seed={} fragment_length={} min_fragment_length={} split_n_run={}",
            manifest.generation_id, manifest.k, manifest.w, manifest.minimizer_hash_seed,
            manifest.fragment_length, manifest.min_fragment_length, manifest.split_n_run)?;
        writeln!(
            writer,
            "reference_id\tshard_index\toriginal_bases\tmapped_bases"
        )?;
        for entry in &manifest.shards {
            let sketch =
                ReferenceSketch::load(&shard_entry_path(prefix, entry), params, false, runtime)?;
            for file in sketch.files {
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}",
                    file.path, entry.shard_index, file.original_length, file.mapped_length
                )?;
            }
        }
        Ok(())
    }
}
