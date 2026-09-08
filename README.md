# fasterANI

`fasterANI` estimates average nucleotide identity (ANI) between query and
reference FASTA genomes. FASTA inputs may be uncompressed or gzip-compressed.

## Install or build

Install the executable from a source checkout:

```bash
cargo install --path .
```

For a repository-local release build instead:

```bash
cargo build --release
```

Executables are written to `target/release/fasterANI` and its lowercase alias
`target/release/fasterani`. Use Rust 1.89 or newer; this change was tested with 1.98.1.

## Quick start

Compare FASTA files directly:

```bash
target/release/fasterANI query \
  --reference references/ref.fna \
  --query queries/query.fna \
  --out results.tsv \
  --threads 8
```

Build an on-disk reference sketch from a file containing one FASTA path per
line. The `sketch` command builds without querying:

```bash
target/release/fasterANI sketch \
  --reference-list references.txt \
  --output sketches/reference-db \
  --threads 8
```

Query that sketch later without resupplying the references:

```bash
target/release/fasterANI query \
  --reference-sketch sketches/reference-db \
  --query queries/query.fna \
  --out results.tsv \
  --threads 8
```

Results are tab-separated and go to stdout unless `--out` is supplied. Runtime
parameters and summaries go to stderr. Run `fasterANI --help` for the complete
CLI and see [the input instructions](docs/input-instructions.md) for the TOML
params-file format. `--threads` limits CPU workers; manifest-backed queries may
also use one bounded background thread to prefetch the next shard.

## Reproducible runs

- Record the ordered reference/query inputs and explicitly set algorithm
  options such as `--kmer-size`, `--window-size`, and
  `--minimizer-hash-seed` (default `42`).
- Keep the startup parameter record from stderr, or store the same values in a
  params file and pass it with `--params-file`.
- Saved-database subcommands inherit sketch settings and reject conflicting overrides.
- Reusing a sketch with reference inputs verifies their ordered path labels and
  count for original builds, or ordered stored identifiers after an update. Neither
  validates FASTA byte contents. Use an explicit remove/add update after editing an
  assembly in place.
- Leave input validation enabled. `--skip-validation` bypasses FASTA path and
  size checks, while `--quiet` suppresses the startup parameter record and final
  summary.

## Update or inspect a saved database

```bash
fasterani inspect --reference-sketch sketches/reference-db
fasterani update --reference-sketch sketches/reference-db --add-list additions.txt
fasterani update --reference-sketch sketches/reference-db --remove-list removals.txt
# Replace assemblies in one transaction:
fasterani update --reference-sketch sketches/reference-db \
  --remove-list old-assemblies.txt --add-list updated-assemblies.txt
```

Addition lists contain FASTA paths (plain or gzip), one per line. Removal lists
contain exact stored reference identifiers, one per line: the FASTA basenames
shown by `inspect`, including extensions. Original FASTAs are not required for
removal. Basenames must be unique; adding an existing identifier requires removing
it in the same transaction. This is identity-based validation, not sequence-based
deduplication or automatic accession/version detection.

Updates retain unaffected shards and repack affected shards from stored minimizers.
Global frequencies are updated by merging counts, and the manifest is published
last. Each old manifest is saved as `<prefix>.<generation>.manifest.json`; its shard,
name-sidecar, and frequency artifacts remain available. Restoring that manifest at
the original `<prefix>.manifest.json` selects the old version; preserve the newer
manifest first if you need it too. No automatic artifact garbage collection is
performed. Readers can finish on the previous generation while an update runs;
cooperating writers are serialized by a filesystem lock.

`sketch` and `update` accept list files only, never positional FASTAs or
`--reference`. `update` can also accept query inputs to query after a successful
update. The previous option-only interface remains available. Updates require a
manifest-backed database; legacy single-file caches remain queryable through the
option-only interface. See [notes.md](notes.md) for the implementation and test summary.
