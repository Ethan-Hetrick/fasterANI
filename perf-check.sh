#!/bin/bash

set -euo pipefail

BINARY=./target/release/fasterANI
SKETCH=~/TESTING/gtdb_genomes_reps_r207-5K-new
GENOME_LIST=~/scratch/gtdb_genomes_reps_r207/genome-list-5K.txt
QUERY=assets/test-data/Escherichia_coli_str_K12_MG1655.fna
TMP=~/scratch/tmp/fasterANI
SHARD_SIZE=22000
RESULTS_DIR=bench_$(date +%Y%m%d_%H%M%S)
REPS=${REPS:-1}

rm -r ~/TESTING/gtdb_genomes_reps_r207-5K-new*

mkdir -p "$RESULTS_DIR"

pick_events() {
    local base_events=(
        major-faults
        minor-faults
        cache-misses
        cache-references
        cycles
        instructions
        branches
        branch-misses
        L1-dcache-loads
        L1-dcache-load-misses
        LLC-loads
        LLC-load-misses
        dTLB-loads
        dTLB-load-misses
        stalled-cycles-frontend
        stalled-cycles-backend
    )

    local optional_events=(
        mem_load_retired.l1_hit
        mem_load_retired.l2_hit
        mem_load_retired.l3_hit
        mem_load_retired.l3_miss
        mem_load_retired.local_dram
        mem_load_retired.remote_dram
    )

    {
        printf "%s\n" "${base_events[@]}"

        for event in "${optional_events[@]}"; do
            if perf list "$event" 2>/dev/null | grep -q "$event"; then
                printf "%s\n" "$event"
            fi
        done
    } | paste -sd, -
}

EVENTS="$(pick_events)"

echo "Results dir: $RESULTS_DIR"
echo "Perf events: $EVENTS"

if [ ! -f "${SKETCH}.manifest.json" ]; then
    echo "Building sketch..."
    "$BINARY" \
        --reference-list "$GENOME_LIST" \
        --reference-sketch "$SKETCH" \
        --shard-size "$SHARD_SIZE" \
        --index-build-mode partitioned \
        --tmp "$TMP" \
        --threads 8 \
        --verbose \
        1>/dev/null
fi

warm_cache() {
    find "$(dirname "$SKETCH")" \
        \( -name "*.fasketch" -o -name "*.fasketch.bgz" \) \
        -type f -print0 |
    xargs -0 cat > /dev/null
}

run_bench() {
    local log=$1
    shift

    /usr/bin/time -v \
    perf stat \
        -d -d -d \
        -e "$EVENTS" \
    "$BINARY" \
        --reference-list "$GENOME_LIST" \
        --query "$QUERY" \
        --reference-sketch "$SKETCH" \
        --shard-size "$SHARD_SIZE" \
        --index-build-mode partitioned \
        --tmp "$TMP" \
        --verbose \
        "$@" \
        1>/dev/null 2>"$log"
}

for THREADS in 1 2 4 8; do
    for REP in $(seq 1 "$REPS"); do
        LOG="$RESULTS_DIR/${THREADS}t_1s_rep${REP}"

        echo "=== threads=$THREADS concurrent_shards=1 rep=$REP ==="

        warm_cache
        run_bench "${LOG}.log" \
            --threads "$THREADS" \
            --max-concurrent-shards 1

        echo "Saved: ${LOG}.log"
    done
done