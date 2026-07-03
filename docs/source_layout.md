# Source layout

The crate is split into a library `src/lib.rs`, crate name `faster_ani`) and a

thin binary `src/main.rs`, which just calls `faster_ani::run()`). All logic

lives in `src/ani/`, organized by concern. Unit tests are co-located with the

code they cover in `#[cfg(test)] mod tests` blocks; end-to-end CLI tests live in

`tests/cli.rs`.

```

src/

  [lib.rs](http://lib.rs)              library crate root; a fast, FastANI-style average nucleotide identity estimator.

  [main.rs](http://main.rs)             binary entry point (faster_ani::run())

  ani/

    [mod.rs](http://mod.rs)            module wiring + re-exports + pub use run

    [constants.rs](http://constants.rs)      algorithm/format constants + type aliases

    [error.rs](http://error.rs)          typed AniError enum (thiserror) + io::Error conversion

    model/

      [reference.rs](http://reference.rs)    reference data types, manifests, ReferenceIndex, param structs

      [query.rs](http://query.rs)        query-side data types: fragments, mapping results, and ANI summaries.

    [runtime.rs](http://runtime.rs)        RuntimeOptions, memory/RSS helpers, progress

    [minimizer.rs](http://minimizer.rs)      canonical minimizer extraction, sliding-window sketches, and query fragmentation.

    [mash.rs](http://mash.rs)           mash-distance math and shared-minimizer lower-bound estimates.

    [metrics.rs](http://metrics.rs)        optional hot-path mapping metrics (debug builds) and seed-hit histograms.

    io_[util.rs](http://util.rs)        path/gzip helpers, byte slicing, ScratchFile

    [mmap.rs](http://mmap.rs)           MmapFile + mmap-backed reference views

    [validation.rs](http://validation.rs)     validation and default-value helpers for CLI and runtime parameters.

    sketch/

      [serialize.rs](http://serialize.rs)    sketch/manifest/sidecar I/O

      [partition.rs](http://partition.rs)    sharding + partition planning

      [build.rs](http://build.rs)        impl ReferenceSketch (in-memory build/index)

      [persist.rs](http://persist.rs)      sketch save/load + streamed cache

      [stream.rs](http://stream.rs)       streaming/partitioned sketch construction

      [lookup.rs](http://lookup.rs)       candidate-region query lookups

      [database.rs](http://database.rs)     SketchDatabase orchestration (collect/load/sharded build)

    [mapping.rs](http://mapping.rs)        seed-hit candidate discovery and sliding-window ANI scoring per query/reference pair.

    [cli.rs](http://cli.rs)            CliArgs + argument parsing

    [pipeline.rs](http://pipeline.rs)       output orchestration + run()

    test_[support.rs](http://support.rs)   shared unit-test fixtures (debug/test only)

tests/

  [cli.rs](http://cli.rs)              end-to-end CLI golden-output test

```