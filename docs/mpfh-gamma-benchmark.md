```bash
set -o pipefail

for gamma in 1.1 3 5 7 13 17 20
do
    run_dir="/run/media/ethan/Lexar/tmp/fasterani-gamma-$gamma"

    echo "Running gamma=$gamma"

    rm -rf -- "$run_dir"
    mkdir -p -- "$run_dir" || continue

    RAYON_NUM_THREADS=4 /usr/bin/time -v -- \
        ./target/performance/fasterANI \
        --reference-list ./benchmarks/ref-1000.txt \
        --reference-sketch "$run_dir/sketch" \
        --mphf-gamma "$gamma" \
        --threads 4 \
        --tmp "$run_dir" \
        --verbose

    status=$?

    if (( status == 0 ))
    then
        printf 'RESULT\tgamma=%s\tsketch_bytes=%s\n' \
            "$gamma" \
            "$(stat -c '%s' "$run_dir/sketch.1.fasketch")"
    else
        printf 'FAILED\tgamma=%s\texit=%s\n' \
            "$gamma" \
            "$status"
    fi

    rm -rf -- "$run_dir"
done 2>&1 | tee gamma-sweep-test.log
```


| Gamma | MPHF build time | Sketch size | Peak RSS |
| ----: | --------------: | ----------: | -------: |
|     3 |         10.41 s |    5.42 GiB | 2.27 GiB |
|     5 |          9.49 s |    5.49 GiB | 2.37 GiB |
|     7 |          9.34 s |    5.55 GiB | 2.48 GiB |
|    10 |      **9.29 s** |    5.64 GiB | 2.60 GiB |
|    13 |          9.42 s |    5.72 GiB | 2.74 GiB |
|    17 |          9.67 s |    5.81 GiB | 2.90 GiB |
|    20 |          9.83 s |    5.88 GiB | 3.02 GiB |
