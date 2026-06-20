# Introduction

`fasterANI` is a reimplementation of [FastANI](https://github.com/ParBLiSS/FastANI) in Rust. `fasterANI` improves upon `fastANI`:

- Supports option to save reference minimizer lookup table for efficient, repeated queries
    - Uses a minimal perfect hash (MPH) to improve hash table size scaling, and subsequently memory scaling and lookup speed
- Addresses ambiguous base handling:
    - Reference minimizers with ambiguous are skipped for candidate scoring
    - Query minimizers with ambiguous bases are no longer factored into fragment identity scoring
    - Supports splitting contigs on a defined number of continuous Ns: `--split-N <int>`
    - The ambiguous base handling allowed for 2-bit encoding of minimizer sequences, which replaced the previous minimizer hash as the lookup key, which slightly improves algorithmic determinism by removing the chance of collisions, thus removing the possibilitity of false positives.
- Full paramaterization of algorithm, with presets `--mode {fast,sensitive,accurate}`, aliases for parameter combinations that improve the use-case applications for the tool.
 - `--mapping-stats` output improves upon the `.visual` output of `fastANI`, giving standard 0-indexed, BED-like format for use in downstream processes, such as re-estimating ANI over coding regions to get `orthoANI`. It contains more fragment-level summary statistics, which can be helpful for detecting horizontal gene transfer, assembly contamination, and genetic hybrids
- Scalability is a hard-requirement of this reimplementation. `fastANI` often looses favor to tools like `skani`, as large queries can explode memory and estimated runtime. The goal of this reimplementation is to scale the more accurate and sensitive `fastANI` algorithm to make it more practical in high-throughput environments
    - Leveraging simd processing significantly speeds up minimizer computations and other algorithm operations, making best use of modern compute resources
    - Certain algorithm modifications improve performance: Skipping over dense ambiguous nucleotide clusters and minimizers, querying a saved reference MPH, optional `--minmer <int>` parameter to use a subset of minimizers for candidate scoring
    - Several options have been implemented to help scale to variable memory, processor and storage settings:
        - `--sketch`: Create reference database once on a high-performance machine (or download one), queries then use minimal resources
        - `--threads`: Standard multi-threading, applies to reference building and querying
        - `--max-memory-gb`: Limits reference shard-size and load-size to prevent peak memory from spiking too high
        - `--index-build-mode <auto|partitioned>` and memory-aware `--shard-minimizers` control reference lookup shard size for tunable memory scaling
        - `--bgzip`: Reads/writes bgzip-compressed saved reference minimizer lookup tables for low storage scenarios
        - `--tmp <path>`: Optionally redirect temporary files to user defined path
        - Multiple algorithmic parameters (e.g. `--kmer-size`, `--kmer-size`, etc.) can be tuned for sensitivity or speed
- Improvements on user-friendliness:
    - Code base is well documented
    - Options like `--header` and `--verbose` improve output and logging readability
    - Syntax is standardized and modular
    - Clear and robust tool documentation through Mkdocs and well-structured CLI `--help` function
    - CI/CD for standardized community contribution and future tool updates
- Adherence, and potential improvement of `fastANI`'s accuracy with respect to `ANIb`
    - Early `fasterANI` testing resulted in certain parameter configurations that resulted in higher pearson correlations to `ANIm`, another gold-standard ANI method, than `fastANI`.
    - Experimental parameter `--min-fragment-length` factors in fragments that are below the designated `--fragment-length <int>`, unlike `fastANI` which exclusively discards these, resulting in poor ANI estimation in low-quality genome pairs (e.g. N50<10k). This may have partially attributed to high ANIm correlations in early test trials.
    - Factoring out ambiguous nucleotides from fragment seeding and identity estimation theoretically aligns better with an ANI-BLAST, which by default doesn't score alignments with ambiguous bases, or treats them as a mismatch
