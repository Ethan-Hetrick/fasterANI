#!/bin/bash
set -euo pipefail

BINARY=./target/release/fasterANI
SKETCH=~/TESTING/gtdb_genomes_reps_r207-fulldb
GENOME_LIST=~/scratch/gtdb_genomes_reps_r207/genome-list.txt
QUERY=assets/test-data/Escherichia_coli_str_K12_MG1655.fna
TMP=~/scratch/tmp/fasterANI
BUILD_THREADS=50
RESULTS_DIR=bench_$(date +%Y%m%d_%H%M%S)
REPS=${REPS:-1}

mkdir -p "$RESULTS_DIR" "$TMP"

rm -r "${SKETCH}"* 2>/dev/null || true

pick_events() {
    local events=(
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
        mem_load_retired.l1_hit
        mem_load_retired.l2_hit
        mem_load_retired.l3_hit
        mem_load_retired.l3_miss
        mem_load_retired.local_dram
        mem_load_retired.remote_dram
    )

    {
        for event in "${events[@]}"; do
            if perf list "$event" 2>/dev/null | grep -q "$event"; then
                printf "%s\n" "$event"
            fi
        done
    } | paste -sd, -
}

EVENTS="$(pick_events)"

echo "Results dir: $RESULTS_DIR"
echo "Perf events: $EVENTS"

run_reference_build() {
    local log="$RESULTS_DIR/reference_build.log"

    echo "=== reference_build threads=$BUILD_THREADS ==="

    /usr/bin/time -v \
    perf stat \
        -d -d -d \
        -e "$EVENTS" \
    "$BINARY" \
        --reference-list "$GENOME_LIST" \
        --reference-sketch "$SKETCH" \
        --index-build-mode partitioned \
        --tmp "$TMP" \
        --threads "$BUILD_THREADS" \
        --verbose \
        1>/dev/null 2>"$log"

    echo "Saved: $log"
}

warm_cache() {
    find "$(dirname "$SKETCH")" \
        \( -name "*.fasketch" \
        -o -name "*.fasketch.bgz" \
        -o -name "*.manifest.json" \
        -o -name "*.contigs.tsv" \
        -o -name "*.contigs.tsv.bgz" \) \
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
        --index-build-mode partitioned \
        --tmp "$TMP" \
        --verbose \
        "$@" \
        1>/dev/null 2>"$log"
}

run_reference_build

for THREADS in 1 8 16 32 50; do
    for REP in $(seq 1 "$REPS"); do
        LOG="$RESULTS_DIR/${THREADS}t_1s_rep${REP}"

        echo "=== query threads=$THREADS concurrent_shards=1 rep=$REP ==="

        warm_cache
        run_bench "${LOG}.log" \
            --threads "$THREADS" \
            --max-concurrent-shards 1

        echo "Saved: ${LOG}.log"
    done
done
