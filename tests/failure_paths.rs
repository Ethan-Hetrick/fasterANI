//! End-to-end checks for failure paths that involve the sharded loader thread.

use std::{
    fs,
    io::{BufReader, Read},
    path::{Path, PathBuf},
    process::{Child, Command, Stdio},
    thread,
    time::{Duration, Instant, SystemTime, UNIX_EPOCH},
};

fn temp_test_dir(name: &str) -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("system clock before unix epoch")
        .as_nanos();
    let path = std::env::temp_dir().join(format!("fasterani-failure-{name}-{nanos}"));
    fs::create_dir_all(&path).expect("create temp test dir");
    path
}

fn write_test_fasta(path: &Path, seed: u64) {
    let mut state = seed;
    let mut sequence = String::with_capacity(6_000);
    for _ in 0..6_000 {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        sequence.push(match state & 3 {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            _ => 'T',
        });
    }
    fs::write(path, format!(">sequence\n{sequence}\n")).expect("write test FASTA");
}

fn terminate_and_collect(mut child: Child, reader: thread::JoinHandle<String>) -> String {
    let _ = child.kill();
    let _ = child.wait();
    reader.join().expect("stderr reader thread panicked")
}

#[test]
fn sharded_loader_error_exits_promptly_without_hanging() {
    let exe = env!("CARGO_BIN_EXE_fasterANI");
    let temp_dir = temp_test_dir("loader-error");
    let reference_one = temp_dir.join("reference-one.fna");
    let reference_two = temp_dir.join("reference-two.fna");
    let sketch_prefix = temp_dir.join("database");
    write_test_fasta(&reference_one, 0x1234_5678_9abc_def0);
    write_test_fasta(&reference_two, 0x0fed_cba9_8765_4321);

    let build_output = Command::new(exe)
        .args([
            "--reference",
            reference_one.to_str().expect("UTF-8 reference path"),
            "--reference",
            reference_two.to_str().expect("UTF-8 reference path"),
            "--reference-sketch",
            sketch_prefix.to_str().expect("UTF-8 sketch prefix"),
            "--max-shard-minimizers",
            "1",
            "--threads",
            "1",
            "--quiet",
        ])
        .output()
        .expect("launch sharded sketch build");
    assert!(
        build_output.status.success(),
        "build failed: {}",
        String::from_utf8_lossy(&build_output.stderr)
    );

    let manifest_path = temp_dir.join("database.manifest.json");
    let manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(&manifest_path).expect("read sharded sketch manifest"))
            .expect("parse sharded sketch manifest");
    let shards = manifest["shards"]
        .as_array()
        .expect("manifest shards array");
    assert_eq!(shards.len(), 2, "test build should create two shards");
    let second_shard = temp_dir.join(
        shards[1]["filename"]
            .as_str()
            .expect("second shard filename"),
    );
    let expected_second_shard_bytes = shards[1]["file_bytes"]
        .as_u64()
        .expect("second shard byte count");
    let mut corrupted_shard = fs::read(&second_shard).expect("read second shard");
    assert_eq!(
        corrupted_shard.len() as u64,
        expected_second_shard_bytes,
        "manifest should record the second shard's exact size"
    );
    corrupted_shard[0] ^= 0xff;
    fs::write(&second_shard, &corrupted_shard).expect("corrupt second shard in place");
    assert_eq!(
        fs::metadata(&second_shard)
            .expect("stat corrupted second shard")
            .len(),
        expected_second_shard_bytes,
        "corruption must preserve size so manifest validation reaches the loader thread"
    );

    let mut child = Command::new(exe)
        .args([
            "--query",
            reference_one.to_str().expect("UTF-8 query path"),
            "--reference-sketch",
            sketch_prefix.to_str().expect("UTF-8 sketch prefix"),
            "--threads",
            "1",
            "--verbose",
            "--quiet",
        ])
        .stdout(Stdio::null())
        .stderr(Stdio::piped())
        .spawn()
        .expect("launch sharded query");
    let stderr = child.stderr.take().expect("child stderr pipe");
    let reader = thread::spawn(move || {
        let mut stderr = BufReader::new(stderr);
        let mut collected = String::new();
        if let Err(error) = stderr.read_to_string(&mut collected) {
            collected.push_str(&format!("stderr read failed: {error}\n"));
        }
        collected
    });

    let deadline = Instant::now() + Duration::from_secs(10);
    let status = loop {
        match child.try_wait().expect("poll sharded query") {
            Some(status) => break status,
            None if Instant::now() < deadline => thread::sleep(Duration::from_millis(25)),
            None => {
                let stderr = terminate_and_collect(child, reader);
                panic!("loader error did not terminate within 10 seconds:\n{stderr}");
            }
        }
    };
    let stderr = reader.join().expect("stderr reader thread panicked");

    assert!(!status.success(), "corrupted shard unexpectedly succeeded");
    assert!(
        stderr.lines().any(|line| {
            line.contains("stage=shard_load")
                && line.contains("event=start")
                && line.contains("shard=2")
        }),
        "query did not reach the second shard in the loader thread: {stderr}"
    );
    assert!(
        stderr.contains("invalid magic header"),
        "corrupted-shard error was not reported: {stderr}"
    );

    let _ = fs::remove_dir_all(temp_dir);
}
