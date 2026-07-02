## Implementation Plan: Pipelined Shard Loading with Parallel Querying

### What needs to change

Only `src/ani/pipeline.rs` requires modification, specifically the `SketchDatabase::Sharded` branch inside the per-query loop. No other files need to change.

### The target architecture

```
loader thread: [load shard 1] → send → [load shard 2] → send → [load shard 3] → ...
main thread:                  ← recv → [query shard 1 on N threads] ← recv → [query shard 2 ...]
```

The channel has a bound of 1, meaning the loader stays at most one shard ahead. At any moment, at most 2 shards are live: one being queried, one being prefetched into the page cache by the loader.

### Step 1: Replace the Rayon parallel shard loop with a `sync_channel`

Remove the existing `pool.install(|| active_shards.par_iter()...)` block entirely. Replace it with a `std::sync::mpsc::sync_channel::<io::Result<ShardQueryResult>>(1)` — the bound of 1 is the key to keeping only one prefetched shard in memory at a time.

The sender side runs on a dedicated thread spawned with `std::thread::spawn`. The receiver side runs on the main thread (or the calling query thread), processing shards one at a time.

### Step 2: The loader thread

```rust
let (tx, rx) = std::sync::mpsc::sync_channel::<io::Result<ShardQueryResult>>(1);

// Clone everything the loader needs before moving into the thread.
let loader_active_shards: Vec<&ShardManifestEntry> = active_shards.clone();
let loader_prefix = prefix.clone();
// ... clone args fields needed: kmer_size, window_size, etc.

std::thread::spawn(move || {
    for (shard_offset, shard) in loader_active_shards.iter().enumerate() {
        // progress logging: event=start

        let load_result = ReferenceSketch::load(
            &shard_entry_path(&loader_prefix, shard),
            SketchParams { ... },
            mapping_stats_requested,
            tmp_dir,
            runtime_options,
        );

        // Send Ok(sketch) or Err — receiver will surface the error.
        // sync_channel blocks here if the receiver hasn't consumed the previous shard yet.
        // That backpressure is exactly the memory bound we want.
        let send_result = tx.send(load_result.map(|sketch| LoadedShard {
            shard_index: shard.shard_index,
            shard_offset,
            sketch,
        }));

        if send_result.is_err() {
            break; // Receiver dropped (error on query side), stop loading.
        }
    }
    // tx drops here, closing the channel.
});
```

Add a small private struct to carry what the receiver needs:

```rust
struct LoadedShard {
    shard_index: usize,
    shard_offset: usize,
    sketch: ReferenceSketch,
}
```

### Step 3: The receiver / query loop

Replace the `for mut shard_result in shard_results` merge loop with a `while let Ok(result) = rx.recv()` loop that both queries and merges inline, avoiding the need to accumulate all shard results before processing:

```rust
while let Ok(load_result) = rx.recv() {
    let loaded: LoadedShard = load_result?; // propagate load errors

    let frequency_threshold = loaded.sketch.index.frequency_threshold(args.freq_threshold_percent);

    // Query using the full thread pool — all args.threads threads available since
    // loading is now off the main thread and not competing.
    let mut raw_stats = collect_query_mappings(
        &loaded.sketch,
        &query_file,
        args.kmer_size,
        args.window_size,
        args.mash_threshold,
        args.mash_confidence,
        args.threads,           // full thread count, not divided
        frequency_threshold,
        performance_metrics_enabled,
    )?;

    let reference_file_offset = reference_file_offsets[loaded.shard_offset];
    let reference_contig_offset = reference_contig_offsets[loaded.shard_offset];

    for mapping in &mut raw_stats.mapping_results {
        mapping.reference_file_id += reference_file_offset;
        mapping.reference_contig_id += reference_contig_offset;
    }

    query_reference_files.append(&mut loaded.sketch.files.into_iter().collect());
    // ... append contig names if mapping_stats_requested
    query_mapping_results.append(&mut raw_stats.mapping_results);
    mapping_elapsed += raw_stats.mapping_elapsed;
    mapping_count += raw_stats.mapping_results.len(); // already moved, use pre-append count

    check_memory_limit("after shard query", runtime_options)?;
    // progress logging: event=complete
}
```

Note that `mapping_threads_per_shard` and `shard_query_parallelism` are no longer needed — the full `args.threads` count is always available to the query since loading is now off the main thread entirely.

### Step 4: Remove the now-unused Rayon shard pool

Delete `shard_query_parallelism`, `mapping_threads_per_shard`, the `AtomicUsize completed_shard_queries`, and the `rayon::ThreadPoolBuilder` that were used for the old parallel-shard approach. The existing query-level Rayon pool inside `map_query_to_reference_parallel` is unaffected and does all the fragment-level parallelism.

### Step 5: `madvise` prefetch hint in `MmapFile::open`

In `src/ani/mmap.rs`, `MmapFile::open` currently issues `MADV_RANDOM`. That's correct for the query phase. The loader thread should issue `MADV_SEQUENTIAL` immediately after `mmap()` and before returning, then the main thread's random access pattern takes over naturally. The simplest approach is to add an optional parameter or a second `pub(crate) fn open_prefetch(path)` variant that calls `MADV_SEQUENTIAL` instead, which the loader uses. Alternatively, add a `prefetch()` method:

```rust
pub(crate) fn prefetch_sequential(&self) {
    unsafe {
        libc::madvise(self.ptr.as_ptr() as *mut libc::c_void, self.len, libc::MADV_SEQUENTIAL);
    }
}
```

Call this on the loaded sketch's mmap immediately after `ReferenceSketch::load()` returns in the loader thread, before sending through the channel. This tells the kernel to start streaming pages into the page cache while the previous shard is still being queried on the main thread.

### Step 6: `max_concurrent_shards` disposition

`max_concurrent_shards` is now meaningless — the channel bound of 1 is the concurrency policy. Leave the CLI flag in place to avoid breaking existing scripts, but ignore its value in the new code path. Add a note in the `--help` text or a startup warning if it's set to a value > 1 indicating that sequential loading is now the policy and this flag has no effect.

### Step 7: Error handling across the thread boundary

The loader sends `io::Result<LoadedShard>`. The receiver propagates errors with `?`. If the query side errors first (e.g. memory limit exceeded), it drops `rx`, which causes the loader's `tx.send()` to return `Err(SendError)`, breaking the loader loop cleanly without a panic. No explicit abort mechanism is needed.

### Summary of changed files

| File | Change |
|------|--------|
| `src/ani/pipeline.rs` | Replace Rayon parallel shard iter with `sync_channel(1)` + loader thread + sequential receiver/query loop; remove `shard_query_parallelism`/`mapping_threads_per_shard`; pass full `args.threads` to `collect_query_mappings` |
| `src/ani/mmap.rs` | Add `prefetch_sequential()` method to `MmapFile`; call it from loader thread after load |
