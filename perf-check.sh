#!/bin/bash

set -euo pipefail

BINARY=./target/release/fasterANI
SKETCH=~/TESTING/gtdb_genomes_reps_r207-5K
GENOME_LIST=~/scratch/gtdb_genomes_reps_r207/genome-list-5K.txt
QUERY=assets/test-data/Escherichia_coli_str_K12_MG1655.fna
TMP=~/scratch/tmp/fasterANI
SHARD_SIZE=22000
RESULTS_DIR=bench_$(date +%Y%m%d_%H%M%S)

mkdir -p "$RESULTS_DIR"

if [ ! -f "${SKETCH}.manifest.json" ]; then
    echo "Building sketch..."
    "$BINARY" \
        --reference-list "$GENOME_LIST" \
        --reference-sketch "$SKETCH" \
        --shard-size $SHARD_SIZE \
        --index-build-mode partitioned \
        --tmp "$TMP" \
        --threads 8 \
        --verbose \
        1>/dev/null
fi

warm_cache() {
    find "$(dirname "$SKETCH")" -name "*.fasketch" -o -name "*.fasketch.bgz" | xargs cat > /dev/null
}

run_bench() {
    local log=$1
    shift
    /usr/bin/time -v \
    perf stat -e major-faults,minor-faults,cache-misses,cache-references \
    "$BINARY" \
        --reference-list "$GENOME_LIST" \
        --query "$QUERY" \
        --reference-sketch "$SKETCH" \
        --shard-size $SHARD_SIZE \
        --index-build-mode partitioned \
        --tmp "$TMP" \
        --verbose \
        "$@" \
        1>/dev/null 2>"$log"
}

for THREADS in 1 2 4 8; do
    for CONCURRENT in 1; do
        LOG="$RESULTS_DIR/${THREADS}t_${CONCURRENT}s"
        echo "=== threads=$THREADS concurrent_shards=$CONCURRENT ==="

        warm_cache
        run_bench "${LOG}.log" --threads "$THREADS" --max-concurrent-shards "$CONCURRENT"
        echo "Saved: ${LOG}.log"
    done
done