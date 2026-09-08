# Input Instructions

`fasterANI` can load runtime options from a TOML file with `--params-file`.
Values in the params file override built-in defaults. Values passed on the CLI
override values from the params file.

Relative paths in the params file are resolved relative to the TOML file's
directory.

## Complete Params File Example

```toml
# Inputs
# Direct FASTA paths. Repeat by adding more entries to the arrays.
reference_files = [
  "references/ref-a.fna",
  "references/ref-b.fna.gz",
]

query_files = [
  "-",
  "queries/query-a.fna",
  "queries/query-b.fna.gz",
]

# Files containing one FASTA path per line. These list files may include blank
# lines and comments beginning with "#"; both are ignored.
reference_lists = [
  "references/reference-list.txt",
]

query_lists = [
  "queries/query-list.txt",
]

# Optional on-disk reference sketch prefix. When queries are omitted, this
# builds the sketch only. When queries are present, the sketch is built or reused.
reference_sketch = "sketches/example-reference-db"

# Label for a streamed stdin query. Use only when query_files contains "-" or
# "/dev/stdin".
query_name = "stdin-query.fna"


# Output
# Omit out to write result TSV to stdout.
out = "results/ani.tsv"

# Prepend a header row to the result TSV.
header = true

# Report ANI independently for each query contig.
per_contig = false

# Write per-fragment mapping statistics.
mapping_stats = "results/mapping-stats.tsv"

# Print progress and diagnostics to stderr.
verbose = false

# Suppress startup parameter summary and final summary logs.
quiet = false


# Seeding
# K-mer size for minimizers. Valid range: 1..=16.
kmer_size = 16

# Minimizer window size. Must be at least 1.
window_size = 24

# Hash seed used for reference and query minimizers. Valid range:
# 0..=4_294_967_295.
minimizer_hash_seed = 42

# Keep only the n smallest-hash minimizers per query fragment for candidate
# scoring. Must be at least 1. Omit this key to use all minimizers.
minmer_count = 1000

# Ignore reference minimizers occurring in more than this percent of reference
# positions. Valid range: 0..=100.
freq_threshold_percent = 0.0


# Fragmenting
# Query fragment length in bp. Must be at least 1.
fragment_length = 3000

# Step between fragment starts. Valid range: 1..=fragment_length.
fragment_stride = 3000

# Keep trailing fragments at least this long. Valid range:
# 1..=fragment_length.
min_fragment_length = 3000

# Split contigs at runs of at least this many ambiguous N bases. Use 0 to
# disable N-run splitting.
split_n_run = 0


# Fragment mapping
# Minimum Mash identity for a query fragment to count toward final ANI.
# Valid range: 0..=100.
mash_threshold = 80.0

# Confidence interval width used by the Mash upper-identity fragment filter.
# Valid range: 0..=1. Use 0 to disable the confidence relaxation; 1 is
# accepted but usually not advised.
mash_confidence = 0.9


# Sketch database / sharding
# Persistent sketch builds use bounded temporary partitions. Omit
# reference_sketch for a one-time, in-memory reference/query run.
# MPHF size/build-time tradeoff for saved sketches. Must be finite and > 1.01.
mphf_gamma = 10.0

# Target maximum persisted size per shard. Accepts raw bytes or B, KiB, MiB,
# and GiB suffixes. A single oversized reference remains in one shard.
max_shard_size = "10GiB"

# Query only selected shard indices. Requires reference_sketch.
# Format matches --shards: comma-separated indices and ranges.
shards = "1,3,5-8"


# Resources
# CPU worker limit. Must be at least 1. Manifest-backed queries may additionally
# use one bounded background thread to prefetch the next shard from storage.
threads = 1

# Directory for temporary shard files.
tmp = "tmp"
```

## Notes

- `skip_validation` is not a params-file option. Use the CLI flag
  `--skip-validation` when you want to bypass runtime FASTA file checks.
- `query_name` labels the stdin query path (`"-"` or `"/dev/stdin"`).
- `minimizer_hash_seed` must match between reference and query minimizers.
- Existing sketches compare supplied references by ordered path label and count;
  use a new sketch prefix after changing a FASTA in place.
- Reference and query list files may include blank lines and comments beginning
  with `#`.
- By default, FASTA paths must exist, be regular files, and be larger than 100
  bytes.
- `minmer_count` is optional; remove it to use all minimizers per query fragment.
- `shards` is optional and only applies when `reference_sketch` is set.

## Startup Parameter Record

Unless `quiet = true` or `--quiet` is used, `fasterANI` writes its complete
effective runtime configuration to stderr as TOML-compatible key/value lines
with source comments. The params-file path is included as commented, CLI-only
metadata; every uncommented key can be copied into a params file. Shard
selections are emitted in the supported string form, for example
`shards = "1,3,5"`. The record always includes validation behavior as commented
CLI-only metadata so provenance stays explicit without adding an unsupported
TOML key.


## Database lifecycle subcommands

Use `query` for one-time comparisons or saved-database searches, `sketch` to build,
`update` to add/remove/replace references, and `inspect` to list stored identifiers.
The lowercase `fasterani` executable is an alias for `fasterANI`. The option-only
CLI remains supported. `--params` is an alias for `--params-file`.

New TOML fields are `command` (query/sketch/update/inspect), `add_lists`, and
`remove_lists`. A CLI subcommand overrides the TOML command. Update list entries
from TOML and CLI are combined. For example:

```toml
command = "update"
reference_sketch = "sketches/reference-db"
add_lists = ["new-assemblies.txt"]
remove_lists = ["old-identifiers.txt"]
threads = 8
```

For a build, use `command = "sketch"`, `reference_lists = ["references.txt"]`,
and `reference_sketch = "sketches/reference-db"`. The CLI `--output` spelling
selects the database prefix for `sketch`; `--out` selects query/inspect TSV output.
Do not use `reference_files` with sketch/update, or `reference_lists` with update.
Each list contains one entry per line; blank lines and lines beginning with `#`
are ignored. Removal identifiers are exact, case-sensitive stored FASTA basenames,
not filesystem paths. Use `inspect` to find them.

Existing path resolution is preserved: CLI reference/add lists resolve relative
FASTA entries from the working directory; lists supplied through TOML resolve
entries relative to the list's directory. The list paths themselves are relative
to the TOML file when specified there. Absolute FASTA paths work in either form.

Saved-database subcommands inherit k, w, minimizer seed, fragment length, minimum
fragment length and N-splitting settings, rejecting conflicting explicit values.
The shard-size target is inherited unless overridden for new additions. Workers,
scratch directory and query-specific settings remain configurable. Startup records
include the command, update list paths and the provenance of inherited parameters.
