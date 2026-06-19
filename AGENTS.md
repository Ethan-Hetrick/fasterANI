# AGENTS.md

Guidance for AI coding agents (and humans) working in this repository.

## What this is

`fasterANI` is a Rust reimplementation of FastANI for average-nucleotide-identity
estimation. It is a **library + binary**: logic lives in the `faster_ani` library
crate (`src/lib.rs` -> `src/ani/`), and `src/main.rs` is a thin wrapper that calls
`faster_ani::run()`.

## Build & test

```bash
cargo build              # debug
cargo build --release    # optimized
cargo test               # unit tests + end-to-end golden CLI test (tests/cli.rs)
cargo fmt --all          # format
cargo clippy --all-targets
./scripts/regression.sh  # full behavior-preservation check (tests + release smoke run)
```

`target-cpu=native` is configured in `.cargo/config.toml`, so you do **not** need
to pass `RUSTFLAGS` manually. (Before that config existed, release builds failed
because the `ensure_simd` dependency hard-errors without AVX2/NEON.) See the
comment in `.cargo/config.toml` for the portability trade-off.

## Gotchas that waste time if you don't know them

- **`cargo` may not be on `PATH` in non-interactive shells.** The Rust toolchain
  lives at `~/.cargo/bin`. Interactive shells get it via `~/.bashrc`, but
  non-interactive shells (e.g. CI steps, some agent shells) source neither
  `~/.bashrc` nor `~/.profile`. If `cargo: command not found`, prepend it:
  `export PATH="$HOME/.cargo/bin:$PATH"`. (RUSTFLAGS is no longer needed thanks to
  `.cargo/config.toml`.)
- **`CARGO_TARGET_DIR` may be redirected.** Some sandboxed environments set
  `CARGO_TARGET_DIR` to a cache directory outside the repo, so freshly built
  artifacts do NOT appear under `./target/`. Do not judge whether a rebuild
  happened by the mtime of `./target/release/fasterANI`. To run the binary you
  just built, prefer `cargo run --release -- <args>` (it always runs the actual
  artifact) instead of invoking `./target/release/fasterANI` by path.
- **Clippy warns about `too_many_arguments`** on the query-mapping / sketch
  functions. These are known, intentional, and non-blocking (warnings, not
  errors). Do not "fix" them by suppressing unless asked.
- **Debug-only code is `#[cfg(debug_assertions)]`-gated** (performance metrics,
  memory accounting). Imports used only by that code must also be gated, or
  release builds emit unused-import warnings. See `pipeline.rs` / `runtime.rs`.

## Conventions

- **No behavior changes without explicit intent.** Output (stdout TSV) must stay
  byte-for-byte stable; verify with `./scripts/regression.sh`. The golden value is
  asserted once in `tests/cli.rs` (single source of truth).
- Unit tests are co-located with the code they cover in `#[cfg(test)] mod tests`
  blocks; shared fixtures live in `src/ani/test_support.rs`.
- Progress (`PROGRESS\t...`) goes to **stderr** and is only printed with
  `--verbose`. Results go to **stdout**.
- A pre-commit hook in `.githooks/pre-commit` runs fmt + clippy + tests. Enable it
  with `git config core.hooksPath .githooks`.

## Source layout

See the "Source layout" section in `README.md` for the per-file breakdown.
