#!/usr/bin/env bash
# Heap profile via DHAT (writes dhat-heap.json in the repo root).
#
# View the JSON at: https://nnethercote.github.io/dh_view/dh_view.html
# (click "Load dhat-heap.json" or drag-and-drop the file).
#
# Uses the `performance` profile so allocations resemble release while keeping symbols.
#
# Usage:
#   ./scripts/dhat.sh --reference ref.fa --query query.fa [extra fasterANI args...]
#   ./scripts/dhat.sh -- --reference ref.fa --query query.fa --threads 1
#
# Default demo (E. coli self-query) when no args are given.
set -euo pipefail

cd "$(dirname "$0")/.."
export PATH="${HOME}/.cargo/bin:${PATH}"

REF="assets/test-data/Escherichia_coli_str_K12_MG1655.fna"
QUERY="assets/test-data/Escherichia_coli_str_K12_MG1655.fna"
OUT="${DHAT_OUT:-dhat-heap.json}"

args=()
if [[ $# -eq 0 ]]; then
    args=(--reference "$REF" --query "$QUERY" --threads 1)
elif [[ "${1:-}" == "--" ]]; then
    shift
    args=("$@")
else
    args=("$@")
fi

rm -f "$OUT"

echo "==> cargo run --profile performance --features dhat-heap"
echo "    args: ${args[*]}"
echo "    output: $OUT"

cargo run --profile performance --features dhat-heap --quiet -- "${args[@]}" >/dev/null

if [[ ! -f "$OUT" ]]; then
    echo "error: $OUT was not created" >&2
    exit 1
fi

bytes="$(wc -c < "$OUT" | tr -d ' ')"
echo "==> wrote $OUT ($bytes bytes)"
echo "    view: https://nnethercote.github.io/dh_view/dh_view.html"
