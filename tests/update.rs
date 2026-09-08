//! Small synthetic end-to-end tests for database lifecycle and CLI policy.
use serde_json::Value;
use std::{
    fs,
    io::Write,
    path::{Path, PathBuf},
    process::{Command, Output},
    time::{SystemTime, UNIX_EPOCH},
};
struct Fixture(PathBuf);
impl Fixture {
    fn new() -> Self {
        let p = std::env::temp_dir().join(format!(
            "fasterani-update-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&p).unwrap();
        Self(p)
    }
    fn path(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }
    fn list(&self, name: &str, entries: &[PathBuf]) -> PathBuf {
        let p = self.path(name);
        fs::write(
            &p,
            entries
                .iter()
                .map(|p| format!("{}\n", p.display()))
                .collect::<String>(),
        )
        .unwrap();
        p
    }
}
impl Drop for Fixture {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}
fn run(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_fasterANI"))
        .args(args)
        .args(["--quiet", "--threads", "2"])
        .output()
        .unwrap()
}
fn ok(args: &[&str]) -> Output {
    let out = run(args);
    assert!(
        out.status.success(),
        "{args:?}\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    out
}
fn bad(args: &[&str], message: &str) {
    let out = run(args);
    assert!(!out.status.success(), "unexpected success: {args:?}");
    assert!(
        String::from_utf8_lossy(&out.stderr).contains(message),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
}
fn txt(p: &Path) -> &str {
    p.to_str().unwrap()
}
fn manifest(db: &Path) -> Value {
    serde_json::from_slice(&fs::read(format!("{}.manifest.json", db.display())).unwrap()).unwrap()
}
fn artifact(db: &Path, name: &str) -> PathBuf {
    db.parent().unwrap().join(name)
}
fn sorted(bytes: Vec<u8>) -> Vec<String> {
    let text = String::from_utf8(bytes).unwrap();
    let mut lines: Vec<_> = text.lines().map(str::to_owned).collect();
    lines.sort();
    lines
}
fn genome(seed: u64, mutate: usize) -> Vec<u8> {
    let mut state = seed;
    let mut bases: Vec<_> = (0..18000)
        .map(|i| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            let base = b"ACGT"[(state & 3) as usize];
            if mutate > 0 && i % mutate == 0 {
                b"ACGT"[((state + 1) & 3) as usize]
            } else {
                base
            }
        })
        .collect();
    bases[5000..5032].fill(b'N');
    let mut out = b">first contig\n".to_vec();
    out.extend_from_slice(&bases[..9000]);
    out.extend_from_slice(b"\n>second contig\n");
    out.extend_from_slice(&bases[9000..]);
    out.push(b'\n');
    out
}
fn build(list: &Path, db: &Path) {
    ok(&[
        "sketch",
        "--reference-list",
        txt(list),
        "--output",
        txt(db),
        "--max-shard-size",
        "1MiB",
        "--kmer-size",
        "13",
        "--window-size",
        "17",
        "--fragment-length",
        "2000",
        "--min-fragment-length",
        "1000",
        "--split-N",
        "20",
    ]);
}
fn compare_fresh(f: &Fixture, db: &Path, references: &[PathBuf], query: &Path, phase: &str) {
    let fresh = f.path(&format!("fresh-{phase}"));
    let list = f.list(&format!("fresh-{phase}.txt"), references);
    build(&list, &fresh);
    for filter in ["0", "1"] {
        let updated_stats = f.path("updated-stats.tsv");
        let fresh_stats = f.path("fresh-stats.tsv");
        let result = ok(&[
            "query",
            "--reference-sketch",
            txt(db),
            "--query",
            txt(query),
            "--max-reference-frequency",
            filter,
            "--mapping-stats",
            txt(&updated_stats),
        ]);
        let expected = ok(&[
            "query",
            "--reference-sketch",
            txt(&fresh),
            "--query",
            txt(query),
            "--max-reference-frequency",
            filter,
            "--mapping-stats",
            txt(&fresh_stats),
        ]);
        assert!(!result.stdout.is_empty(), "empty comparison in {phase}");
        assert_eq!(
            sorted(result.stdout),
            sorted(expected.stdout),
            "ANI mismatch: {phase} filter={filter}"
        );
        assert_eq!(
            sorted(fs::read(updated_stats).unwrap()),
            sorted(fs::read(fresh_stats).unwrap()),
            "mapping mismatch: {phase}"
        );
    }
    let a = manifest(db);
    let b = manifest(&fresh);
    assert_eq!(
        fs::read(artifact(
            db,
            a["global_frequency_filename"].as_str().unwrap()
        ))
        .unwrap(),
        fs::read(artifact(
            &fresh,
            b["global_frequency_filename"].as_str().unwrap()
        ))
        .unwrap(),
        "frequency mismatch: {phase}"
    );
}
#[test]
fn updates_match_fresh_builds_and_preserve_old_generations() {
    let f = Fixture::new();
    let db = f.path("db");
    let a = f.path("a.fa");
    let b = f.path("b.fa");
    let c = f.path("c.fa.gz");
    let query = f.path("query.fa");
    fs::write(&a, genome(42, 0)).unwrap();
    fs::write(&b, genome(42, 113)).unwrap();
    fs::write(&query, genome(42, 0)).unwrap();
    let mut gz = flate2::write::GzEncoder::new(
        fs::File::create(&c).unwrap(),
        flate2::Compression::default(),
    );
    gz.write_all(&genome(42, 97)).unwrap();
    gz.finish().unwrap();
    let initial = f.list("initial.txt", &[a.clone(), b.clone()]);
    build(&initial, &db);
    let mut old = manifest(&db);
    // The upstream perf manifest predates the optional identifier inventory.
    old.as_object_mut().unwrap().remove("reference_identifiers");
    fs::write(
        format!("{}.manifest.json", db.display()),
        serde_json::to_vec_pretty(&old).unwrap(),
    )
    .unwrap();
    assert_eq!(old["shards"].as_array().unwrap().len(), 1);
    let old_path = artifact(&db, old["shards"][0]["filename"].as_str().unwrap());
    let old_bytes = fs::read(&old_path).unwrap();
    let old_manifest_bytes = fs::read(format!("{}.manifest.json", db.display())).unwrap();
    let additions = f.list("add.txt", std::slice::from_ref(&c));
    ok(&[
        "update",
        "--reference-sketch",
        txt(&db),
        "--add-list",
        txt(&additions),
    ]);
    let added = manifest(&db);
    assert_eq!(added["total_references"], 3);
    assert_eq!(added["shards"].as_array().unwrap().len(), 2);
    assert_eq!(fs::read(&old_path).unwrap(), old_bytes);
    assert_eq!(
        fs::read(format!(
            "{}.{}.manifest.json",
            db.display(),
            old["generation_id"].as_str().unwrap()
        ))
        .unwrap(),
        old_manifest_bytes
    );
    compare_fresh(&f, &db, &[a.clone(), b.clone(), c.clone()], &query, "add");
    let c_shard = artifact(&db, added["shards"][1]["filename"].as_str().unwrap());
    let c_bytes = fs::read(&c_shard).unwrap();
    // Fail after an affected shard has been repacked but before publication.
    let bad_gzip = f.path("bad.fa.gz");
    fs::write(&bad_gzip, vec![0xff; 1000]).unwrap();
    let bad_add = f.list("bad-add.txt", &[bad_gzip]);
    let removal = f.path("remove.txt");
    fs::write(&removal, "a.fa\n").unwrap();
    let current_bytes = fs::read(format!("{}.manifest.json", db.display())).unwrap();
    let failure = run(&[
        "update",
        "--reference-sketch",
        txt(&db),
        "--remove-list",
        txt(&removal),
        "--add-list",
        txt(&bad_add),
    ]);
    assert!(!failure.status.success());
    assert_eq!(
        fs::read(format!("{}.manifest.json", db.display())).unwrap(),
        current_bytes
    );
    assert_eq!(fs::read(&old_path).unwrap(), old_bytes);
    // Original FASTAs may be gone. A copy is kept only for the independent fresh-build oracle.
    fs::create_dir(f.path("oracle")).unwrap();
    let saved_b = f.path("oracle/b.fa");
    fs::copy(&b, &saved_b).unwrap();
    fs::remove_file(&a).unwrap();
    fs::remove_file(&b).unwrap();
    ok(&[
        "update",
        "--reference-sketch",
        txt(&db),
        "--remove-list",
        txt(&removal),
    ]);
    assert_eq!(fs::read(&c_shard).unwrap(), c_bytes);
    assert_eq!(fs::read(&old_path).unwrap(), old_bytes);
    compare_fresh(&f, &db, &[saved_b, c.clone()], &query, "remove");
    fs::create_dir(f.path("v2")).unwrap();
    let new_b = f.path("v2/b.fa");
    fs::write(&new_b, genome(42, 71)).unwrap();
    let replacements = f.list("replace.txt", std::slice::from_ref(&new_b));
    fs::write(&removal, "b.fa\n").unwrap();
    // Same stored identifier is permitted only with explicit removal in this transaction.
    bad(
        &[
            "update",
            "--reference-sketch",
            txt(&db),
            "--add-list",
            txt(&replacements),
        ],
        "already exists",
    );
    let replacement_params = f.path("replacement.toml");
    fs::write(&replacement_params,"command = \"update\"\nreference_sketch = \"db\"\nadd_lists = [\"replace.txt\"]\nremove_lists = [\"remove.txt\"]\n").unwrap();
    let immediate = ok(&["--params", txt(&replacement_params), "--query", txt(&query)]);
    assert!(!immediate.stdout.is_empty());
    assert_eq!(fs::read(&c_shard).unwrap(), c_bytes);
    compare_fresh(&f, &db, &[c.clone(), new_b.clone()], &query, "replace");
    let inspected = ok(&["inspect", "--reference-sketch", txt(&db)]);
    assert!(String::from_utf8_lossy(&inspected.stdout).contains("b.fa\t"));
    let remove_c = f.path("remove-c.txt");
    fs::write(&remove_c, "c.fa.gz\n").unwrap();
    let params = f.path("update.toml");
    fs::write(
        &params,
        "command = \"update\"\nreference_sketch = \"db\"\nremove_lists = [\"remove-c.txt\"]\n",
    )
    .unwrap();
    ok(&["--params-file", txt(&params)]);
    compare_fresh(
        &f,
        &db,
        std::slice::from_ref(&new_b),
        &query,
        "params-remove",
    );
    // Empty database is a valid update result, and can subsequently be repopulated.
    fs::write(&removal, "b.fa\n").unwrap();
    ok(&[
        "update",
        "--reference-sketch",
        txt(&db),
        "--remove-list",
        txt(&removal),
    ]);
    assert_eq!(manifest(&db)["total_references"], 0);
    let empty = ok(&[
        "query",
        "--reference-sketch",
        txt(&db),
        "--query",
        txt(&query),
    ]);
    assert!(empty.stdout.is_empty());
    ok(&[
        "update",
        "--reference-sketch",
        txt(&db),
        "--add-list",
        txt(&replacements),
    ]);
    compare_fresh(&f, &db, &[new_b], &query, "repopulate");
    // Restore a snapshot at the same prefix: all original artifacts remain usable.
    fs::write(
        format!("{}.manifest.json", db.display()),
        &old_manifest_bytes,
    )
    .unwrap();
    let restored = ok(&[
        "query",
        "--reference-sketch",
        txt(&db),
        "--query",
        txt(&query),
    ]);
    let restored = String::from_utf8(restored.stdout).unwrap();
    assert!(restored.contains("a.fa"));
    assert!(restored.contains("b.fa"));
}
#[test]
fn subcommands_enforce_list_inputs_parameters_and_writer_lock() {
    let f = Fixture::new();
    let a = f.path("a.fa");
    fs::write(&a, genome(9, 0)).unwrap();
    let db = f.path("db");
    let list = f.list("refs.txt", std::slice::from_ref(&a));
    bad(
        &["sketch", "--reference", txt(&a), "--output", txt(&db)],
        "not supported",
    );
    bad(
        &["sketch", txt(&a), "--output", txt(&db)],
        "unknown argument",
    );
    bad(
        &[
            "update",
            "--reference",
            txt(&a),
            "--reference-sketch",
            txt(&db),
        ],
        "not supported",
    );
    bad(
        &[
            "update",
            "--reference-sketch",
            txt(&db),
            "--add-list",
            txt(&list),
        ],
        "existing",
    );
    build(&list, &db);
    bad(
        &[
            "sketch",
            "--reference-list",
            txt(&list),
            "--output",
            txt(&db),
        ],
        "already exists",
    );
    bad(
        &[
            "query",
            "--reference-sketch",
            txt(&db),
            "--query",
            txt(&a),
            "--kmer-size",
            "16",
        ],
        "conflicting",
    );
    bad(&["update", "--reference-sketch", txt(&db)], "nonempty");
    let removal = f.path("remove.txt");
    fs::write(&removal, "unknown.fa\n").unwrap();
    bad(
        &[
            "update",
            "--reference-sketch",
            txt(&db),
            "--remove-list",
            txt(&removal),
        ],
        "unknown removal",
    );
    fs::write(&removal, "a.fa\na.fa\n").unwrap();
    bad(
        &[
            "update",
            "--reference-sketch",
            txt(&db),
            "--remove-list",
            txt(&removal),
        ],
        "duplicate removal",
    );
    fs::write(&removal, "a.fa\n").unwrap();
    let lock = fs::OpenOptions::new()
        .read(true)
        .write(true)
        .open(format!("{}.lock", db.display()))
        .unwrap();
    lock.lock().unwrap();
    bad(
        &[
            "update",
            "--reference-sketch",
            txt(&db),
            "--remove-list",
            txt(&removal),
        ],
        "writer lock",
    );
    drop(lock);
    let p = f.path("invalid.toml");
    fs::write(
        &p,
        format!("command = \"sketch\"\nreference_files = [{:?}]\n", txt(&a)),
    )
    .unwrap();
    bad(&["--params-file", txt(&p)], "list files");
    let one = ok(&["query", "--reference-list", txt(&list), "--query", txt(&a)]);
    assert!(!one.stdout.is_empty());
}
