#!/usr/bin/env bash
# CPU flamegraph via `cargo flamegraph` + Linux perf.
#
# Uses the `performance` profile (release opts + debug symbols) from Cargo.toml.
# Output: flamegraph.svg in the repo root (override with FLAMEGRAPH_OUT=path.svg).
#
# Requires: cargo-flamegraph (`cargo install flamegraph`), perf, write access to perf events.
# On WSL/Linux, if perf is blocked, lower paranoid temporarily (needs root):
#   sudo sysctl -w kernel.perf_event_paranoid=1
#
# Usage:
#   ./scripts/flamegraph.sh --reference ref.fa --query query.fa [extra fasterANI args...]
#   ./scripts/flamegraph.sh -- --reference ref.fa --query query.fa --threads 1
#
# Default demo (E. coli self-query) when no args are given.
set -euo pipefail

cd "$(dirname "$0")/.."
export PATH="${HOME}/.cargo/bin:${PATH}"

REF="assets/test-data/Escherichia_coli_str_K12_MG1655.fna"
QUERY="assets/test-data/Escherichia_coli_str_K12_MG1655.fna"
OUT="${FLAMEGRAPH_OUT:-flamegraph.svg}"

if ! command -v cargo-flamegraph >/dev/null 2>&1; then
    echo "error: cargo-flamegraph not found; run: cargo install flamegraph" >&2
    exit 1
fi

if ! command -v perf >/dev/null 2>&1; then
    echo "error: perf not found" >&2
    echo "  Ubuntu/WSL: sudo apt install linux-tools-common linux-tools-\$(uname -r)" >&2
    echo "  then retry; if perf still missing, try: sudo apt install linux-tools-generic" >&2
    exit 1
fi

paranoid="$(cat /proc/sys/kernel/perf_event_paranoid 2>/dev/null || echo unknown)"
if [[ "$paranoid" != unknown && "$paranoid" -gt 1 ]]; then
    echo "warning: kernel.perf_event_paranoid=$paranoid (often need <= 1 for user perf)" >&2
    echo "         try: sudo sysctl -w kernel.perf_event_paranoid=1" >&2
fi

args=()
if [[ $# -eq 0 ]]; then
    args=(--reference "$REF" --query "$QUERY" --threads 1)
elif [[ "${1:-}" == "--" ]]; then
    shift
    args=("$@")
else
    args=("$@")
fi

echo "==> cargo flamegraph --profile performance --output $OUT"
echo "    args: ${args[*]}"

cargo flamegraph --profile performance --output "$OUT" -- "${args[@]}"

echo "==> wrote $OUT (open in browser or https://firefox-devtools.github.io/profiler/)"
