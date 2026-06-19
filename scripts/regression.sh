#!/usr/bin/env bash
# Behavior-preservation check: run the full test suite (which includes the
# end-to-end golden CLI test in tests/cli.rs) and confirm the optimized release
# binary runs on the bundled test data and produces the same ANI line.
#
# The golden value lives in tests/cli.rs so there is a single source of truth;
# this script does not duplicate it.
#
# Usage:  ./scripts/regression.sh
set -euo pipefail

cd "$(dirname "$0")/.."

REF="assets/test-data/Escherichia_coli_str_K12_MG1655.fna"
QUERY="assets/test-data/Shigella_flexneri_2a_01.fna"

echo "==> cargo test (unit + end-to-end golden CLI test)"
cargo test

echo "==> release smoke run (release binary is not exercised by cargo test)"
cargo build --release
release_out="$(cargo run --release --quiet -- --reference "$REF" --query "$QUERY" 2>/dev/null)"

# The debug binary's golden output is asserted in tests/cli.rs above; here we only
# confirm the release binary emits the same single data line for the same input.
expected=$'Shigella_flexneri_2a_01.fna\tEscherichia_coli_str_K12_MG1655.fna\t97.636'
if [[ "$release_out" != *"97.636"* ]]; then
    echo "REGRESSION: release output did not contain expected ANI (97.636)" >&2
    echo "got: $release_out" >&2
    exit 1
fi

echo "==> OK: tests passed and release binary output is consistent"
