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

The executable is written to `target/release/fasterANI`.

## Quick start

Compare FASTA files directly:

```bash
target/release/fasterANI \
  --reference references/ref.fna \
  --query queries/query.fna \
  --out results.tsv \
  --threads 8
```

Build an on-disk reference sketch from a file containing one FASTA path per
line. Omitting query inputs selects build-only mode:

```bash
target/release/fasterANI \
  --reference-list references.txt \
  --reference-sketch sketches/reference-db \
  --threads 8
```

Query that sketch later without resupplying the references:

```bash
target/release/fasterANI \
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
- Use the same seeding settings when building and querying a reference sketch.
- Reusing a sketch with reference inputs verifies their ordered path labels and
  count, not the FASTA byte contents. Rebuild with `--force` after editing a
  reference in place.
- Leave input validation enabled. `--skip-validation` bypasses FASTA path and
  size checks, while `--quiet` suppresses the startup parameter record and final
  summary.
