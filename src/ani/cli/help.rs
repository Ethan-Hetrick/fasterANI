//! Command-line help text.

pub(super) fn usage() -> &'static str {
    "usage: fasterANI [(--reference <ref.fa> | --reference-list <refs.txt>)...] \
[(--query <query.fa> | --query-list <queries.txt>)...] [options]

Commands:
  fasterANI query --reference-list refs.txt --query-list queries.txt
  fasterANI query --reference-sketch db --query-list queries.txt
  fasterANI sketch --reference-list refs.txt --output db
  fasterANI update --reference-sketch db --add-list additions.txt --remove-list removals.txt
  fasterANI inspect --reference-sketch db [--out references.tsv]

  sketch accepts reference lists only; update accepts add/remove lists only.
  Each list contains one entry per line. Removal entries are exact stored FASTA
  basenames (see inspect); the old FASTA files need not exist. Either update list
  may be omitted. Both lists together perform an atomic replacement operation.
  Existing shards are retained; only shards affected by removal are rebuilt.
  Saved database parameters are inherited; conflicting overrides are rejected.
  update optionally accepts queries to run after successful publication.
  The older option-only invocation remains supported.

Inputs:
  --params-file, --params <path>         Load runtime parameters from a TOML file.
                                 Scalar CLI values override scalar file values.
                                 Reference/query inputs from both sources are combined.
                                 Relative paths in the file are resolved relative
                                 to the TOML file's directory.
  --reference <path>           Reference FASTA (optionally gzip-compressed).
  --reference-list <path>      File of reference FASTA paths, one per line.
  --reference-sketch <prefix>  Build/reuse an on-disk reference sketch at this prefix.
                                 Query inputs optional (build-only when omitted).
                                 Persistent builds use bounded temporary partitions.
  --query <path>               Query FASTA (repeatable). Use `-` to read one query from
                                 stdin (optionally gzip-compressed).
  --query-list <path>          File of query FASTA paths, one per line.
  --query-name <label>         Display/path label for a stdin query (with `--query -`).
  --skip-validation            When specified, file checks are not performed on input FASTA files.

Output:
  --out <path>                 Write results TSV here (default: stdout).
  --header                     Prepend a column-name header row to the results TSV
                                 default: off.
  --per-contig                 Report one row per query contig and reference file.
                                 Aggregate genome-pair summaries are written first
                                 as commented lines.
  --mapping-stats <path>       Write a per-fragment mapping-stats TSV (always headered).
  --verbose                    Print PROGRESS/diagnostics to stderr (default: off).
  --quiet, --silent            Suppress startup parameter summary and final SUMMARY
                                 log (default: off).

  Results columns (tab-separated):
    query_file           Query genome file path.
    reference_file       Reference genome file path.
    ANI                  Average nucleotide identity (%).
    AF                   Aligned fraction of query fragments.
    total_fragments      Mappable query bases / fragment_length (non-integer).
    median_ANI           Median fragment ANI; less sensitive to outliers than the ANI.
    stddev               Standard deviation of fragment ANI.
    MAD                  Median absolute deviation of fragment ANI from median_ANI.
    ci_95_upper          95% upper confidence interval for ANI
    ci_95_lower          95% lower confidence interval for ANI
    F99                  Fraction of retained fragments with ANI >= 99%.
    F80                  Fraction of retained fragments with ANI <= 80%.

  Note: Fragment counts may be fractional.
        With --per-contig, contigs that have no usable fragments are reported
        with NaN ANI fields.
        --per-contig output columns are: query_file, reference_file,
        query_contig, eligible_fragments, shared_fragments, shared_bases,
        ANI, median_ANI, stddev, MAD, ci_95_upper, ci_95_lower, F99, F80.

Seeding (minimizer sketch; applies to both references and queries):
  --kmer-size <1..=16>         K-mer size for minimizers (default 16).
  --window-size <n>            Minimizer window size, >= 1 (default 24).
  --minimizer-hash-seed <0..=4294967295>
                               Hash seed for minimizers (default 42).
  --minmer-count <n>           Keep only the n smallest-hash minimizers ('minmers') per query
                                 fragment for candidate scoring; n >= 1
                                 (default behavior uses all).
  --max-reference-frequency <0..=100>
                               Ignore reference minimizers occurring in more than this
                                 percent of reference positions; filters out frequent,
                                 uninformative k-mers (default 0).

Fragmenting (how each query contig is cut into fragments):
  --fragment-length <bp>       Query fragment length, >= 1 (default 3000).
  --fragment-stride <bp>       Step between fragment starts; 1..=fragment-length
                                 default: equal to fragment-length, i.e. non-overlapping.
  --min-fragment-length <bp>   Keep trailing fragments at least this long; 1..=fragment-length
                                 alias: --min-fraglen; default: fragment-length.
  --split-N <bp>               Split contigs at runs of >= this many ambiguous (N)
                               bases; 0 disables splitting (alias: --split-n; default 0).

Fragment mapping (thresholds applied to each individual fragment alignment):
  --mash-threshold <0..=100>   Drop query fragments whose estimated Mash identity is below
                                 this percent before calculating ANI (default 80).
                                 Use 0 to disable this fragment identity filter.
  --mash-confidence <0..=1>    Confidence interval width for the Mash upper-identity bound
                                 used in fragment filtering (default 0.9).
                                 Higher values are more permissive; 0.9 uses one-sided
                                 tail alpha 0.05. Use 0 to require the estimated identity
                                 itself to pass --mash-threshold; 1 is accepted but
                                 usually not advised.

  Per genome pair, fasterANI keeps only reciprocal-best fragment mappings and reports
    ANI as the length-weighted mean of those retained fragments' identities.

Sketch database / sharding:
  --max-shard-size <size>       Target maximum persisted size per shard (default 10GiB).
                                 Accepts raw bytes or B, KiB, MiB, and GiB suffixes.
  --shards <list>              Comma-separated shard indices and ranges to query,
                                 e.g. 1,3,5-8. Queries all shards when omitted.
                                 Requires --reference-sketch.
  --mphf-gamma <float>         MPHF size/build-time tradeoff for saved sketches
                                 (must be finite and > 1.01; default 10).

Resources:
  --threads <n>                CPU worker limit, >= 1 (default 1). Sharded queries may also use
                                 one bounded I/O-prefetch thread.
  --tmp <dir>                  Directory for temporary shard files.
  -h, --help                   Show help.
  -v, --version                Show version."
}
