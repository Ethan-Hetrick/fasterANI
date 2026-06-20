# Installation



```bash
# Build
export RUSTFLAGS="-C target-cpu=native"
cargo build --release

# Help menu
./target/release/fasterANI --help

# Run test data
./target/release/fasterANI 
    --reference assets/test-data/Escherichia_coli_str_K12_MG1655.fna \
    --query assets/test-data/Shigella_flexneri_2a_01.fna 2> /dev/null
```

## Source layout

The crate is split into a library (`src/lib.rs`, crate name `faster_ani`) and a
thin binary (`src/main.rs`, which just calls `faster_ani::run()`). All logic
lives in `src/ani/`, organized by concern. Unit tests are co-located with the
code they cover in `#[cfg(test)] mod tests` blocks; end-to-end CLI tests live in
`tests/cli.rs`.

```
src/
  lib.rs              library crate root (pub mod ani; pub use ani::run)
  main.rs             binary entry point (faster_ani::run())
  ani/
    mod.rs            module wiring + re-exports + pub use run
    constants.rs      algorithm/format constants + type aliases
    error.rs          typed AniError enum (thiserror) + io::Error conversion
    model/
      reference.rs    reference data types, manifests, ReferenceIndex, param structs
      query.rs        query/mapping data types
    runtime.rs        RuntimeOptions, memory/RSS helpers, progress
    minimizer.rs      sliding sketch, Fenwick, canonical minimizers, fragmenting
    mash.rs           mash-distance + shared-minimizer estimates
    metrics.rs        MappingMetrics (debug/release variants), histograms
    io_util.rs        path/gzip helpers, byte slicing, ScratchFile
    mmap.rs           MmapFile + mmap-backed reference views
    validation.rs     validate_*/default_* parameter helpers
    sketch/
      serialize.rs    sketch/manifest/sidecar I/O
      partition.rs    sharding + partition planning
      build.rs        impl ReferenceSketch (in-memory build/index)
      persist.rs      sketch save/load + streamed cache
      stream.rs       streaming/partitioned sketch construction
      lookup.rs       candidate-region query lookups
      database.rs     SketchDatabase orchestration (collect/load/sharded build)
    mapping.rs        seed-hit mapping engine + ANI computation
    cli.rs            CliArgs + argument parsing
    pipeline.rs       output orchestration + run()
    test_support.rs   shared unit-test fixtures (debug/test only)
tests/
  cli.rs              end-to-end CLI golden-output test
```
