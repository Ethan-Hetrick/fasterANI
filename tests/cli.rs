//! End-to-end test that runs the compiled `fasterANI` binary against the bundled
//! test genomes and checks the emitted ANI line. This guards the public CLI
//! contract (arguments in, TSV out) the same way the README example does.

use std::{
    fs,
    io::Write,
    path::{Path, PathBuf},
    process::{Command, Stdio},
    time::{SystemTime, UNIX_EPOCH},
};

fn temp_test_dir(name: &str) -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("system clock before unix epoch")
        .as_nanos();
    let path = std::env::temp_dir().join(format!("fasterani-cli-{name}-{nanos}"));
    fs::create_dir_all(&path).expect("create temp test dir");
    path
}

fn fixture_path(name: &str) -> String {
    std::env::current_dir()
        .expect("current dir")
        .join("assets/test-data")
        .join(name)
        .to_string_lossy()
        .into_owned()
}

fn write_synthetic_fasta(path: &Path, name: &str, seed: u64, len: usize) {
    let mut state = seed;
    let sequence: String = (0..len)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            char::from(b"ACGT"[(state & 3) as usize])
        })
        .collect();
    fs::write(path, format!(">{name}\n{sequence}\n")).expect("write synthetic FASTA");
}

fn assert_expected_test_data_result(stdout: &str) {
    let fields: Vec<&str> = stdout.trim_end().split('\t').collect();
    assert_eq!(fields.len(), 12, "unexpected result fields: {fields:?}");
    assert_eq!(fields[0], "assets/test-data/Shigella_flexneri_2a_01.fna");
    assert_eq!(fields[1], "Escherichia_coli_str_K12_MG1655.fna");
    assert_eq!(
        &fields[2..10],
        &["97.636", "0.807", "1608.00", "98.318", "2.500", "0.645", "97.772", "97.500"]
    );

    for (index, field) in fields[10..].iter().enumerate() {
        field
            .parse::<f64>()
            .unwrap_or_else(|err| panic!("field {} was not numeric: {field:?}: {err}", index + 10));
    }
}

#[test]
fn cli_emits_expected_ani_for_test_data() {
    let exe = env!("CARGO_BIN_EXE_fasterANI");
    let output = Command::new(exe)
        .args([
            "--reference",
            "assets/test-data/Escherichia_coli_str_K12_MG1655.fna",
            "--query",
            "assets/test-data/Shigella_flexneri_2a_01.fna",
        ])
        .output()
        .expect("failed to launch fasterANI binary");

    assert!(
        output.status.success(),
        "binary exited with status {:?}",
        output.status
    );

    let stdout = String::from_utf8(output.stdout).expect("stdout was not valid UTF-8");
    assert_expected_test_data_result(&stdout);
}

#[test]
fn cli_accepts_streamed_query_from_stdin() {
    let exe = env!("CARGO_BIN_EXE_fasterANI");
    let query_bytes =
        fs::read("assets/test-data/Shigella_flexneri_2a_01.fna").expect("read query fixture");

    let mut child = Command::new(exe)
        .args([
            "--reference",
            "assets/test-data/Escherichia_coli_str_K12_MG1655.fna",
            "--query",
            "-",
            "--query-name",
            "assets/test-data/Shigella_flexneri_2a_01.fna",
        ])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("failed to launch fasterANI binary");

    child
        .stdin
        .take()
        .expect("stdin pipe")
        .write_all(&query_bytes)
        .expect("write query FASTA to stdin");

    let output = child
        .wait_with_output()
        .expect("failed to wait for fasterANI binary");

    assert!(
        output.status.success(),
        "binary exited with status {:?}: {}",
        output.status,
        String::from_utf8_lossy(&output.stderr)
    );

    let stdout = String::from_utf8(output.stdout).expect("stdout was not valid UTF-8");
    assert_expected_test_data_result(&stdout);
}

#[test]
fn cli_query_name_overrides_params_and_is_independent_of_argument_order() {
    let exe = env!("CARGO_BIN_EXE_fasterANI");
    let temp_dir = temp_test_dir("query-name-precedence");
    let params_path = temp_dir.join("params.toml");
    let reference = fixture_path("Escherichia_coli_str_K12_MG1655.fna");
    let query = fixture_path("Shigella_flexneri_2a_01.fna");
    let query_bytes = fs::read(&query).expect("read query fixture");
    fs::write(
        &params_path,
        format!(
            r#"
reference_files = ["{reference}"]
query_files = ["-"]
query_name = "from-params"
"#,
        ),
    )
    .expect("write params file");

    let mut child = Command::new(exe)
        .args([
            "--params-file",
            params_path.to_str().expect("utf-8 params path"),
            "--query-name",
            "from-cli",
            "--query",
            &query,
            "--quiet",
        ])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("failed to launch fasterANI binary");

    child
        .stdin
        .take()
        .expect("stdin pipe")
        .write_all(&query_bytes)
        .expect("write query FASTA to stdin");
    let output = child
        .wait_with_output()
        .expect("failed to wait for fasterANI binary");

    assert!(
        output.status.success(),
        "binary exited with status {:?}: {}",
        output.status,
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8(output.stdout).expect("stdout was not valid UTF-8");
    assert!(
        stdout.lines().any(|line| line.starts_with("from-cli\t")),
        "CLI query label was not applied: {stdout}"
    );
    assert!(!stdout.contains("from-params\t"));

    let _ = fs::remove_dir_all(temp_dir);
}

#[test]
fn no_arguments_displays_help() {
    let exe = env!("CARGO_BIN_EXE_fasterANI");
    let output = Command::new(exe)
        .output()
        .expect("failed to launch fasterANI binary");

    assert!(output.status.success());
    let stderr = String::from_utf8(output.stderr).expect("stderr was not valid UTF-8");
    assert!(stderr.starts_with("usage: fasterANI"));
    assert!(stderr.contains("--reference <path>"));
}

#[test]
fn help_documents_hash_seed_and_build_only_diagnostic_names_real_flag() {
    let exe = env!("CARGO_BIN_EXE_fasterANI");
    let help_output = Command::new(exe)
        .arg("--help")
        .output()
        .expect("failed to launch fasterANI binary");
    assert!(help_output.status.success());
    let help = String::from_utf8(help_output.stderr).expect("stderr was not valid UTF-8");
    assert!(help.contains("--minimizer-hash-seed <0..=4294967295>"));
    assert!(help.contains("Hash seed for minimizers (default 42)"));

    let error_output = Command::new(exe)
        .args([
            "--reference",
            "assets/test-data/Escherichia_coli_str_K12_MG1655.fna",
        ])
        .output()
        .expect("failed to launch fasterANI binary");
    assert!(!error_output.status.success());
    let error = String::from_utf8(error_output.stderr).expect("stderr was not valid UTF-8");
    assert!(error.contains("omit queries when using --reference-sketch for build-only mode"));
    assert!(!error.contains("using --sketch for build-only mode"));
}

#[test]
fn shards_flag_requires_reference_sketch() {
    let exe = env!("CARGO_BIN_EXE_fasterANI");
    let output = Command::new(exe)
        .args([
            "--reference",
            "assets/test-data/Escherichia_coli_str_K12_MG1655.fna",
            "--query",
            "assets/test-data/Shigella_flexneri_2a_01.fna",
            "--shards",
            "1",
        ])
        .output()
        .expect("failed to launch fasterANI binary");

    assert!(
        !output.status.success(),
        "binary unexpectedly succeeded: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stderr = String::from_utf8(output.stderr).expect("stderr was not valid UTF-8");
    assert!(stderr.contains("--shards requires --reference-sketch"));
}

#[test]
fn max_reference_frequency_flag_is_accepted() {
    let exe = env!("CARGO_BIN_EXE_fasterANI");
    let temp_dir = temp_test_dir("max-reference-frequency");
    let out_path = temp_dir.join("out.tsv");
    let output = Command::new(exe)
        .args([
            "--reference",
            "assets/test-data/Escherichia_coli_str_K12_MG1655.fna",
            "--query",
            "assets/test-data/Shigella_flexneri_2a_01.fna",
            "--max-reference-frequency",
            "1.5",
            "--out",
            out_path.to_str().expect("utf-8 out path"),
        ])
        .output()
        .expect("failed to launch fasterANI binary");

    assert!(
        output.status.success(),
        "binary exited with status {:?}: {}",
        output.status,
        String::from_utf8_lossy(&output.stderr)
    );
    let stderr = String::from_utf8(output.stderr).expect("stderr was not valid UTF-8");
    assert!(stderr.contains("freq_threshold_percent = 1.5  # from CLI"));

    let _ = fs::remove_dir_all(temp_dir);
}

#[test]
fn direct_saved_and_sharded_workflows_are_exactly_reproducible() {
    let exe = env!("CARGO_BIN_EXE_fasterANI");
    let temp_dir = temp_test_dir("workflow-parity");
    let reference_a_path = temp_dir.join("reference-a.fna");
    let reference_b_path = temp_dir.join("reference-b.fna");
    let query_a_path = temp_dir.join("query-a.fna");
    let query_b_path = temp_dir.join("query-b.fna");

    let mut state: u64 = 0x9e37_79b9_7f4a_7c15;
    let mut reference_a = Vec::with_capacity(12_000);
    for _ in 0..12_000 {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        reference_a.push(b"ACGT"[(state & 3) as usize]);
    }
    let mut reference_b = reference_a.clone();
    for base in reference_b.iter_mut().skip(17).step_by(37) {
        *base = match *base {
            b'A' => b'C',
            b'C' => b'G',
            b'G' => b'T',
            _ => b'A',
        };
    }
    let mut query_a = reference_a.clone();
    for base in query_a.iter_mut().skip(29).step_by(113) {
        *base = match *base {
            b'A' => b'G',
            b'C' => b'T',
            b'G' => b'A',
            _ => b'C',
        };
    }
    let mut query_b = reference_b.clone();
    for base in query_b.iter_mut().skip(41).step_by(97) {
        *base = match *base {
            b'A' => b'T',
            b'C' => b'A',
            b'G' => b'C',
            _ => b'G',
        };
    }

    fs::write(
        &reference_a_path,
        format!(
            ">reference-a\n{}\n",
            String::from_utf8(reference_a).unwrap()
        ),
    )
    .expect("write reference A");
    fs::write(
        &reference_b_path,
        format!(
            ">reference-b\n{}\n",
            String::from_utf8(reference_b).unwrap()
        ),
    )
    .expect("write reference B");
    fs::write(
        &query_a_path,
        format!(">query-a\n{}\n", String::from_utf8(query_a).unwrap()),
    )
    .expect("write query A");
    fs::write(
        &query_b_path,
        format!(">query-b\n{}\n", String::from_utf8(query_b).unwrap()),
    )
    .expect("write query B");

    let reference_a = reference_a_path.to_string_lossy().into_owned();
    let reference_b = reference_b_path.to_string_lossy().into_owned();
    let query_a = query_a_path.to_string_lossy().into_owned();
    let query_b = query_b_path.to_string_lossy().into_owned();
    let common_args = vec![
        "--reference".to_owned(),
        reference_a,
        "--reference".to_owned(),
        reference_b,
        "--query".to_owned(),
        query_a,
        "--query".to_owned(),
        query_b,
        "--kmer-size".to_owned(),
        "8".to_owned(),
        "--window-size".to_owned(),
        "12".to_owned(),
        "--fragment-length".to_owned(),
        "1000".to_owned(),
        "--max-reference-frequency".to_owned(),
        "1".to_owned(),
        "--quiet".to_owned(),
    ];
    let run = |extra_args: &[&str]| {
        let mut args = common_args.clone();
        args.extend(extra_args.iter().map(|arg| (*arg).to_owned()));
        let output = Command::new(exe)
            .args(&args)
            .output()
            .expect("failed to launch fasterANI binary");
        assert!(
            output.status.success(),
            "binary exited with status {:?} for {args:?}: {}",
            output.status,
            String::from_utf8_lossy(&output.stderr)
        );
        String::from_utf8(output.stdout).expect("stdout was not valid UTF-8")
    };
    let read_persisted_artifacts = |prefix: &Path| {
        let manifest_path = PathBuf::from(format!("{}.manifest.json", prefix.display()));
        let manifest: serde_json::Value =
            serde_json::from_slice(&fs::read(manifest_path).expect("read shard manifest"))
                .expect("parse shard manifest");
        let parent: &Path = prefix.parent().expect("artifact directory");
        let shards: Vec<Vec<u8>> = manifest["shards"]
            .as_array()
            .expect("manifest shards array")
            .iter()
            .map(|shard| {
                let filename = shard["filename"].as_str().expect("shard filename");
                fs::read(parent.join(filename)).expect("read persisted shard")
            })
            .collect();
        let frequency_filename = manifest["global_frequency_filename"]
            .as_str()
            .expect("global frequency filename");
        let frequency =
            fs::read(parent.join(frequency_filename)).expect("read global frequency artifact");
        let sidecar_paths: Vec<PathBuf> = fs::read_dir(parent)
            .expect("read artifact directory")
            .map(|entry| entry.expect("artifact directory entry").path())
            .filter(|path| {
                path.file_name()
                    .and_then(|name| name.to_str())
                    .is_some_and(|name| name.starts_with("fasterani-names."))
            })
            .collect();
        assert_eq!(
            sidecar_paths.len(),
            1,
            "expected one content-addressed sidecar"
        );
        let sidecar = fs::read(&sidecar_paths[0]).expect("read contig sidecar");
        (shards, frequency, sidecar)
    };

    let direct_threads_one = run(&["--threads", "1"]);
    assert!(
        !direct_threads_one.is_empty(),
        "synthetic fixture produced no ANI results"
    );
    let query_labels: Vec<&str> = direct_threads_one
        .lines()
        .map(|line| line.split('\t').next().expect("query label column"))
        .collect();
    let first_query_b: usize = query_labels
        .iter()
        .position(|label| label.ends_with("query-b.fna"))
        .expect("query B results");
    assert!(
        first_query_b > 0
            && query_labels[..first_query_b]
                .iter()
                .all(|label| label.ends_with("query-a.fna"))
            && query_labels[first_query_b..]
                .iter()
                .all(|label| label.ends_with("query-b.fna")),
        "query result blocks were missing or out of input order: {query_labels:?}"
    );
    let direct_threads_two = run(&["--threads", "2"]);
    assert_eq!(direct_threads_one, direct_threads_two);

    let mapping_stats_one = temp_dir.join("mapping-stats-threads-1.tsv");
    let mapping_stats_two = temp_dir.join("mapping-stats-threads-2.tsv");
    assert_eq!(
        direct_threads_one,
        run(&[
            "--threads",
            "1",
            "--mapping-stats",
            mapping_stats_one
                .to_str()
                .expect("utf-8 mapping stats path"),
        ])
    );
    assert_eq!(
        direct_threads_one,
        run(&[
            "--threads",
            "2",
            "--mapping-stats",
            mapping_stats_two
                .to_str()
                .expect("utf-8 mapping stats path"),
        ])
    );
    assert_eq!(
        fs::read(&mapping_stats_one).expect("read serial mapping stats"),
        fs::read(&mapping_stats_two).expect("read parallel mapping stats")
    );

    let hash_prefix_one = temp_dir.join("hash-threads-1/database");
    let hash_output_one = run(&[
        "--threads",
        "1",
        "--reference-sketch",
        hash_prefix_one.to_str().expect("utf-8 hash prefix"),
        "--index-build-mode",
        "hash",
    ]);
    assert_eq!(direct_threads_one, hash_output_one);

    let hash_prefix = temp_dir.join("hash-threads-2/database");
    let hash_output = run(&[
        "--threads",
        "2",
        "--reference-sketch",
        hash_prefix.to_str().expect("utf-8 hash prefix"),
        "--index-build-mode",
        "hash",
    ]);
    assert_eq!(direct_threads_one, hash_output);
    assert_eq!(
        read_persisted_artifacts(&hash_prefix_one),
        read_persisted_artifacts(&hash_prefix),
        "serial and parallel hash builds produced different persisted artifacts"
    );

    let partitioned_prefix_one = temp_dir.join("partitioned-threads-1/database");
    let partitioned_output_one = run(&[
        "--threads",
        "1",
        "--reference-sketch",
        partitioned_prefix_one
            .to_str()
            .expect("utf-8 partitioned prefix"),
        "--index-build-mode",
        "partitioned",
    ]);
    assert_eq!(direct_threads_one, partitioned_output_one);

    let partitioned_prefix = temp_dir.join("partitioned-threads-2/database");
    let partitioned_output = run(&[
        "--threads",
        "2",
        "--reference-sketch",
        partitioned_prefix
            .to_str()
            .expect("utf-8 partitioned prefix"),
        "--index-build-mode",
        "partitioned",
    ]);
    assert_eq!(direct_threads_one, partitioned_output);
    assert_eq!(
        read_persisted_artifacts(&partitioned_prefix_one),
        read_persisted_artifacts(&partitioned_prefix),
        "serial and parallel partitioned builds produced different persisted artifacts"
    );

    let sharded_prefix = temp_dir.join("saved-sharded");
    let sharded_output = run(&[
        "--threads",
        "2",
        "--reference-sketch",
        sharded_prefix.to_str().expect("utf-8 sharded prefix"),
        "--index-build-mode",
        "hash",
        "--max-shard-minimizers",
        "1",
    ]);
    assert_eq!(direct_threads_one, sharded_output);

    let manifest_path = PathBuf::from(format!("{}.manifest.json", sharded_prefix.display()));
    let manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(manifest_path).expect("read shard manifest"))
            .expect("parse shard manifest");
    assert_eq!(
        manifest["shards"]
            .as_array()
            .expect("manifest shards array")
            .len(),
        2,
        "max-shard-minimizers=1 should place the two references in separate shards"
    );

    let _ = fs::remove_dir_all(temp_dir);
}

#[test]
fn skipped_queries_warn_count_and_advance_progress_monotonically() {
    let exe = env!("CARGO_BIN_EXE_fasterANI");
    let temp_dir = temp_test_dir("skipped-query-progress");
    let reference_path = temp_dir.join("reference.fna");
    let query_path = temp_dir.join("query.fna");
    let skipped_path = temp_dir.join("all-ambiguous.fna");
    let mut state: u64 = 0xd1b5_4a32_d192_ed03;
    let sequence: String = (0..1_200)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            char::from(b"ACGT"[(state & 3) as usize])
        })
        .collect();
    fs::write(&reference_path, format!(">reference\n{sequence}\n")).expect("write reference");
    fs::write(&query_path, format!(">query\n{sequence}\n")).expect("write query");
    fs::write(&skipped_path, format!(">skipped\n{}\n", "N".repeat(1_200)))
        .expect("write skipped query");

    let output = Command::new(exe)
        .args([
            "--reference",
            reference_path.to_str().expect("utf-8 reference path"),
            "--query",
            query_path.to_str().expect("utf-8 query path"),
            "--query",
            skipped_path.to_str().expect("utf-8 skipped query path"),
            "--kmer-size",
            "8",
            "--window-size",
            "12",
            "--fragment-length",
            "300",
            "--threads",
            "2",
            "--verbose",
        ])
        .output()
        .expect("launch fasterANI");
    assert!(
        output.status.success(),
        "binary exited with status {:?}: {}",
        output.status,
        String::from_utf8_lossy(&output.stderr)
    );
    let stderr = String::from_utf8(output.stderr).expect("stderr UTF-8");
    assert!(stderr.contains("WARNING\tevent=query_skipped\treason=no_usable_fragments"));
    assert!(
        stderr.contains("queries_requested=2\tqueries_processed=1\tqueries_skipped=1"),
        "missing skipped-query summary counts: {stderr}"
    );
    let completion_counts: Vec<usize> = stderr
        .lines()
        .filter(|line| {
            line.contains("stage=query")
                && (line.contains("event=skipped") || line.contains("event=complete"))
        })
        .map(|line| {
            line.split('\t')
                .find_map(|field| field.strip_prefix("query_done="))
                .expect("query_done field")
                .parse::<usize>()
                .expect("numeric query_done")
        })
        .collect();
    assert_eq!(completion_counts, vec![1, 2]);

    let _ = fs::remove_dir_all(temp_dir);
}

#[test]
fn old_freq_threshold_percent_flag_is_rejected() {
    let exe = env!("CARGO_BIN_EXE_fasterANI");
    let output = Command::new(exe)
        .args([
            "--reference",
            "assets/test-data/Escherichia_coli_str_K12_MG1655.fna",
            "--query",
            "assets/test-data/Shigella_flexneri_2a_01.fna",
            "--freq-threshold-percent",
            "1.5",
        ])
        .output()
        .expect("failed to launch fasterANI binary");

    assert!(
        !output.status.success(),
        "binary unexpectedly succeeded: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stderr = String::from_utf8(output.stderr).expect("stderr was not valid UTF-8");
    assert!(stderr.contains("unknown argument \"--freq-threshold-percent\""));
}

#[test]
fn version_short_flag_is_dash_v_not_bare_v() {
    let exe = env!("CARGO_BIN_EXE_fasterANI");
    let version_output = Command::new(exe)
        .arg("-v")
        .output()
        .expect("failed to launch fasterANI binary");
    assert!(
        version_output.status.success(),
        "-v unexpectedly failed: {}",
        String::from_utf8_lossy(&version_output.stderr)
    );
    assert!(String::from_utf8(version_output.stderr)
        .expect("stderr was not valid UTF-8")
        .contains("fasterANI "));

    let bare_v_output = Command::new(exe)
        .arg("v")
        .output()
        .expect("failed to launch fasterANI binary");
    assert!(
        !bare_v_output.status.success(),
        "bare v unexpectedly succeeded"
    );
    assert!(String::from_utf8(bare_v_output.stderr)
        .expect("stderr was not valid UTF-8")
        .contains("unknown argument \"v\""));
}

#[test]
fn reference_is_optional_when_querying_existing_sharded_sketch() {
    let exe = env!("CARGO_BIN_EXE_fasterANI");
    let temp_dir = temp_test_dir("query-existing-sharded-sketch");
    let sketch_prefix = temp_dir.join("database");
    let reference = "assets/test-data/Escherichia_coli_str_K12_MG1655.fna";
    let query = "assets/test-data/Shigella_flexneri_2a_01.fna";

    let build_output = Command::new(exe)
        .args([
            "--reference",
            reference,
            "--reference-sketch",
            sketch_prefix.to_str().expect("utf-8 sketch prefix"),
            "--quiet",
        ])
        .output()
        .expect("failed to launch fasterANI binary");
    assert!(
        build_output.status.success(),
        "build exited with status {:?}: {}",
        build_output.status,
        String::from_utf8_lossy(&build_output.stderr)
    );

    let query_output = Command::new(exe)
        .args([
            "--query",
            query,
            "--reference-sketch",
            sketch_prefix.to_str().expect("utf-8 sketch prefix"),
            "--shards",
            "1",
            "--quiet",
        ])
        .output()
        .expect("failed to launch fasterANI binary");
    assert!(
        query_output.status.success(),
        "query exited with status {:?}: {}",
        query_output.status,
        String::from_utf8_lossy(&query_output.stderr)
    );
    let stderr = String::from_utf8(query_output.stderr).expect("stderr was not valid UTF-8");
    assert!(!stderr.contains("missing --reference"));

    let _ = fs::remove_dir_all(temp_dir);
}

#[test]
fn build_only_existing_reference_sketch_requires_force() {
    let exe = env!("CARGO_BIN_EXE_fasterANI");
    let temp_dir = temp_test_dir("existing-sketch-build-only-requires-force");
    let sketch_prefix = temp_dir.join("database");
    let reference_list = temp_dir.join("references.txt");
    let reference = fixture_path("Escherichia_coli_str_K12_MG1655.fna");
    fs::write(&reference_list, format!("{reference}\n")).expect("write reference list");

    let build_output = Command::new(exe)
        .args([
            "--reference-list",
            reference_list.to_str().expect("utf-8 reference list"),
            "--reference-sketch",
            sketch_prefix.to_str().expect("utf-8 sketch prefix"),
            "--quiet",
        ])
        .output()
        .expect("failed to launch fasterANI binary");
    assert!(
        build_output.status.success(),
        "build exited with status {:?}: {}",
        build_output.status,
        String::from_utf8_lossy(&build_output.stderr)
    );

    let blocked_output = Command::new(exe)
        .args([
            "--reference-list",
            reference_list.to_str().expect("utf-8 reference list"),
            "--reference-sketch",
            sketch_prefix.to_str().expect("utf-8 sketch prefix"),
            "--quiet",
        ])
        .output()
        .expect("failed to launch fasterANI binary");
    assert!(
        !blocked_output.status.success(),
        "build unexpectedly succeeded: {}",
        String::from_utf8_lossy(&blocked_output.stderr)
    );
    let stderr = String::from_utf8(blocked_output.stderr).expect("stderr was not valid UTF-8");
    assert!(stderr.contains(
        "Reference sketch exists and no query was provided. Use `--force` to overwrite. Exiting..."
    ));

    let forced_output = Command::new(exe)
        .args([
            "--reference-list",
            reference_list.to_str().expect("utf-8 reference list"),
            "--reference-sketch",
            sketch_prefix.to_str().expect("utf-8 sketch prefix"),
            "--force",
            "--quiet",
        ])
        .output()
        .expect("failed to launch fasterANI binary");
    assert!(
        forced_output.status.success(),
        "forced build exited with status {:?}: {}",
        forced_output.status,
        String::from_utf8_lossy(&forced_output.stderr)
    );

    let _ = fs::remove_dir_all(temp_dir);
}

#[test]
fn sharded_query_loads_each_reference_shard_once_for_multiple_queries() {
    let exe = env!("CARGO_BIN_EXE_fasterANI");
    let temp_dir = temp_test_dir("multi-query-shard-loads-once");
    let sketch_prefix = temp_dir.join("database");
    let mapping_stats_path = temp_dir.join("mapping-stats.tsv");
    let reference_one_path = temp_dir.join("reference-one.fna");
    let reference_two_path = temp_dir.join("reference-two.fna");
    write_synthetic_fasta(&reference_one_path, "reference-one", 0x1234_5678, 6_000);
    write_synthetic_fasta(&reference_two_path, "reference-two", 0x8765_4321, 6_000);
    let sketch_prefix = sketch_prefix.to_str().expect("utf-8 sketch prefix");
    let mapping_stats_path = mapping_stats_path
        .to_str()
        .expect("utf-8 mapping stats path");
    let reference_one = reference_one_path.to_str().expect("utf-8 reference path");
    let reference_two = reference_two_path.to_str().expect("utf-8 reference path");

    let build_output = Command::new(exe)
        .args([
            "--reference",
            reference_one,
            "--reference",
            reference_two,
            "--reference-sketch",
            sketch_prefix,
            "--max-shard-minimizers",
            "1",
            "--verbose",
            "--quiet",
        ])
        .output()
        .expect("failed to launch fasterANI binary");
    assert!(
        build_output.status.success(),
        "build exited with status {:?}: {}",
        build_output.status,
        String::from_utf8_lossy(&build_output.stderr)
    );
    let build_stderr =
        String::from_utf8(build_output.stderr).expect("build stderr was not valid UTF-8");
    let nested_build_events: Vec<&str> = build_stderr
        .lines()
        .filter(|line| {
            line.contains("stage=reference_build")
                || line.contains("stage=sketch_save")
                || line.contains("stage=sketch_load")
        })
        .collect();
    assert!(!nested_build_events.is_empty());
    for event in nested_build_events {
        assert!(
            event.contains("generation_id="),
            "unattributed event: {event}"
        );
        assert!(event.contains("shard="), "unattributed event: {event}");
        assert!(
            event.contains("assigned_threads="),
            "unattributed event: {event}"
        );
    }

    let query_output = Command::new(exe)
        .args([
            "--query",
            reference_one,
            "--query",
            reference_two,
            "--reference-sketch",
            sketch_prefix,
            "--mapping-stats",
            mapping_stats_path,
            "--verbose",
            "--quiet",
        ])
        .output()
        .expect("failed to launch fasterANI binary");
    assert!(
        query_output.status.success(),
        "query exited with status {:?}: {}",
        query_output.status,
        String::from_utf8_lossy(&query_output.stderr)
    );

    let stdout = String::from_utf8(query_output.stdout).expect("stdout was not valid UTF-8");
    assert!(!stdout.is_empty());
    let stderr = String::from_utf8(query_output.stderr).expect("stderr was not valid UTF-8");
    assert_eq!(stderr.matches("stage=shard_load\tevent=start").count(), 2);
    assert_eq!(
        stderr.matches("stage=shard_load\tevent=complete").count(),
        2
    );
    assert!(!stderr.contains("stage=shard_load\tevent=start\tquery_done="));

    let mapping_stats = fs::read_to_string(mapping_stats_path).expect("read mapping stats");
    assert!(mapping_stats.starts_with("query_file\treference_file\tquery_contig"));
    assert!(mapping_stats.lines().count() > 1);

    let _ = fs::remove_dir_all(temp_dir);
}

#[test]
fn quiet_suppresses_startup_summary_from_cli() {
    let exe = env!("CARGO_BIN_EXE_fasterANI");
    let output = Command::new(exe)
        .args([
            "--reference",
            "assets/test-data/Escherichia_coli_str_K12_MG1655.fna",
            "--query",
            "assets/test-data/Shigella_flexneri_2a_01.fna",
            "--threads",
            "2",
            "--quiet",
        ])
        .output()
        .expect("failed to launch fasterANI binary");

    assert!(
        output.status.success(),
        "binary exited with status {:?}: {}",
        output.status,
        String::from_utf8_lossy(&output.stderr)
    );

    let stdout = String::from_utf8(output.stdout).expect("stdout was not valid UTF-8");
    assert_expected_test_data_result(&stdout);

    let stderr = String::from_utf8(output.stderr).expect("stderr was not valid UTF-8");
    assert!(!stderr.contains("FasterANI effective runtime parameters"));
    assert!(!stderr.contains("threads = 2"));
}

#[test]
fn quiet_suppresses_startup_summary_from_params_file() {
    let exe = env!("CARGO_BIN_EXE_fasterANI");
    let temp_dir = temp_test_dir("params-quiet");
    let params_path = temp_dir.join("params.toml");
    let reference = fixture_path("Escherichia_coli_str_K12_MG1655.fna");
    let query = fixture_path("Shigella_flexneri_2a_01.fna");
    fs::write(
        &params_path,
        format!(
            r#"
reference_files = ["{reference}"]
query_files = ["{query}"]
threads = 4
quiet = true
"#,
        ),
    )
    .expect("write params file");

    let output = Command::new(exe)
        .args([
            "--params-file",
            params_path.to_str().expect("utf-8 params path"),
        ])
        .output()
        .expect("failed to launch fasterANI binary");

    assert!(
        output.status.success(),
        "binary exited with status {:?}: {}",
        output.status,
        String::from_utf8_lossy(&output.stderr)
    );

    let stdout = String::from_utf8(output.stdout).expect("stdout was not valid UTF-8");
    assert!(stdout.contains("\t97.636\t0.807\t1608.00\t"));

    let stderr = String::from_utf8(output.stderr).expect("stderr was not valid UTF-8");
    assert!(!stderr.contains("FasterANI effective runtime parameters"));
    assert!(!stderr.contains("threads = 4"));

    let _ = fs::remove_dir_all(temp_dir);
}

#[test]
fn params_file_provides_defaults() {
    let exe = env!("CARGO_BIN_EXE_fasterANI");
    let temp_dir = temp_test_dir("params-defaults");
    let params_path = temp_dir.join("params.toml");
    let out_path = temp_dir.join("out.tsv");
    let reference = fixture_path("Escherichia_coli_str_K12_MG1655.fna");
    let query = fixture_path("Shigella_flexneri_2a_01.fna");
    fs::write(
        &params_path,
        format!(
            r#"
reference_files = ["{reference}"]
query_files = ["{query}"]
threads = 4
minimizer_hash_seed = 7
per_contig = true
"#,
        ),
    )
    .expect("write params file");

    let output = Command::new(exe)
        .args([
            "--params-file",
            params_path.to_str().expect("utf-8 params path"),
            "--out",
            out_path.to_str().expect("utf-8 out path"),
        ])
        .output()
        .expect("failed to launch fasterANI binary");

    assert!(
        output.status.success(),
        "binary exited with status {:?}: {}",
        output.status,
        String::from_utf8_lossy(&output.stderr)
    );
    let stderr = String::from_utf8(output.stderr).expect("stderr was not valid UTF-8");
    assert!(stderr.contains("threads = 4  # from params file"));
    assert!(stderr.contains("minimizer_hash_seed = 7  # from params file"));
    assert!(stderr.contains("per_contig = true  # from params file"));

    let _ = fs::remove_dir_all(temp_dir);
}

#[test]
fn startup_parameter_record_is_supported_toml() {
    let exe = env!("CARGO_BIN_EXE_fasterANI");
    let temp_dir = temp_test_dir("startup-record");
    let params_path = temp_dir.join("params.toml");
    let missing_sketch = temp_dir.join("missing-sketch");
    let query = fixture_path("Shigella_flexneri_2a_01.fna");
    fs::write(
        &params_path,
        format!(
            r#"
reference_sketch = "{}"
query_files = ["{query}"]
minimizer_hash_seed = 7
shards = "3,1"
force = true
"#,
            missing_sketch.display(),
        ),
    )
    .expect("write params file");

    let output = Command::new(exe)
        .args([
            "--params-file",
            params_path.to_str().expect("utf-8 params path"),
        ])
        .output()
        .expect("failed to launch fasterANI binary");
    let stderr = String::from_utf8(output.stderr).expect("stderr was not valid UTF-8");
    let header = "################## FasterANI effective runtime parameters ##################";
    let record_start = stderr
        .find(header)
        .unwrap_or_else(|| panic!("startup record missing: {stderr}"));
    let record_body = &stderr[record_start + header.len()..];
    let record_body = record_body
        .strip_prefix('\n')
        .expect("startup record header should end with a newline");
    let record_end = record_body
        .find("############################################################################")
        .expect("startup record footer");
    let record = &record_body[..record_end];
    let parsed: toml::Value = toml::from_str(record).expect("startup record should be valid TOML");
    let table = parsed.as_table().expect("startup record should be a table");

    assert!(!table.contains_key("params_file"));
    assert_eq!(table["minimizer_hash_seed"].as_integer(), Some(7));
    assert_eq!(table["shards"].as_str(), Some("1,3"));
    assert_eq!(table["force"].as_bool(), Some(true));
    assert_eq!(table["threads"].as_integer(), Some(1));
    assert_eq!(table["per_contig"].as_bool(), Some(false));
    assert_eq!(table["index_build_mode"].as_str(), Some("auto"));
    for required_key in [
        "header",
        "verbose",
        "quiet",
        "freq_threshold_percent",
        "kmer_size",
        "window_size",
        "fragment_length",
        "fragment_stride",
        "min_fragment_length",
        "mash_threshold",
        "mash_confidence",
        "mphf_gamma",
        "split_n_run",
        "max_shard_minimizers",
    ] {
        assert!(
            table.contains_key(required_key),
            "effective configuration omitted {required_key}: {record}"
        );
    }
    assert!(record.contains("# params_file = "));
    assert!(record.contains("# skip_validation = false  (CLI-only metadata)"));
    assert!(record.contains("minimizer_hash_seed = 7  # from params file"));
    assert!(record.contains("shards = \"1,3\"  # from params file"));

    let replay_path = temp_dir.join("replay.toml");
    fs::write(&replay_path, record).expect("write replay params file");
    let replay_output = Command::new(exe)
        .args([
            "--params-file",
            replay_path.to_str().expect("utf-8 replay path"),
        ])
        .output()
        .expect("failed to replay startup record");
    let replay_stderr =
        String::from_utf8(replay_output.stderr).expect("stderr was not valid UTF-8");
    assert!(
        !replay_stderr.contains("Failed to parse params file"),
        "replayed startup record did not parse: {replay_stderr}"
    );

    let _ = fs::remove_dir_all(temp_dir);
}

#[test]
fn cli_args_override_params_file() {
    let exe = env!("CARGO_BIN_EXE_fasterANI");
    let temp_dir = temp_test_dir("params-override");
    let params_path = temp_dir.join("params.toml");
    let out_path = temp_dir.join("out.tsv");
    let reference = fixture_path("Escherichia_coli_str_K12_MG1655.fna");
    let query = fixture_path("Shigella_flexneri_2a_01.fna");
    fs::write(
        &params_path,
        format!(
            r#"
reference_files = ["{reference}"]
query_files = ["{query}"]
threads = 4
minimizer_hash_seed = 7
"#,
        ),
    )
    .expect("write params file");

    let output = Command::new(exe)
        .args([
            "--params-file",
            params_path.to_str().expect("utf-8 params path"),
            "--threads",
            "2",
            "--minimizer-hash-seed",
            "8",
            "--out",
            out_path.to_str().expect("utf-8 out path"),
        ])
        .output()
        .expect("failed to launch fasterANI binary");

    assert!(
        output.status.success(),
        "binary exited with status {:?}: {}",
        output.status,
        String::from_utf8_lossy(&output.stderr)
    );
    let stderr = String::from_utf8(output.stderr).expect("stderr was not valid UTF-8");
    assert!(stderr.contains("threads = 2  # from CLI"));
    assert!(!stderr.contains("threads = 4  # from params file"));
    assert!(stderr.contains("minimizer_hash_seed = 8  # from CLI"));
    assert!(!stderr.contains("minimizer_hash_seed = 7  # from params file"));

    let _ = fs::remove_dir_all(temp_dir);
}

#[test]
fn params_file_and_cli_references_combine() {
    let exe = env!("CARGO_BIN_EXE_fasterANI");
    let temp_dir = temp_test_dir("params-combine");
    let params_path = temp_dir.join("params.toml");
    let out_path = temp_dir.join("out.tsv");
    let ref_from_file = temp_dir.join("ref-from-file.fna");
    let ref_from_cli = temp_dir.join("ref-from-cli.fna");
    let query = temp_dir.join("query.fna");
    let sequence = "ACGT".repeat(1000);
    fs::write(&ref_from_file, format!(">ref_from_file\n{sequence}\n")).expect("write ref");
    fs::write(&ref_from_cli, format!(">ref_from_cli\n{sequence}\n")).expect("write ref");
    fs::write(&query, format!(">query\n{sequence}\n")).expect("write query");
    fs::write(
        &params_path,
        format!(
            r#"
reference_files = ["{}"]
query_files = ["{}"]
"#,
            ref_from_file.display(),
            query.display(),
        ),
    )
    .expect("write params file");

    let output = Command::new(exe)
        .args([
            "--params-file",
            params_path.to_str().expect("utf-8 params path"),
            "--reference",
            ref_from_cli.to_str().expect("utf-8 ref path"),
            "--out",
            out_path.to_str().expect("utf-8 out path"),
        ])
        .output()
        .expect("failed to launch fasterANI binary");

    let stderr = String::from_utf8(output.stderr).expect("stderr was not valid UTF-8");
    assert!(stderr.contains("reference_files = ["));
    assert!(stderr.contains("ref-from-file.fna"));
    assert!(stderr.contains("ref-from-cli.fna"));
    assert!(stderr.contains("# from params file + CLI"));

    let _ = fs::remove_dir_all(temp_dir);
}

#[test]
fn params_file_references_are_validated() {
    let exe = env!("CARGO_BIN_EXE_fasterANI");
    let temp_dir = temp_test_dir("params-invalid-reference");
    let params_path = temp_dir.join("params.toml");
    let query = fixture_path("Shigella_flexneri_2a_01.fna");
    fs::write(
        &params_path,
        format!(
            r#"
reference_files = ["missing-reference.fna"]
query_files = ["{query}"]
"#,
        ),
    )
    .expect("write params file");

    let output = Command::new(exe)
        .args([
            "--params-file",
            params_path.to_str().expect("utf-8 params path"),
        ])
        .output()
        .expect("failed to launch fasterANI binary");

    assert!(
        !output.status.success(),
        "binary unexpectedly succeeded: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stderr = String::from_utf8(output.stderr).expect("stderr was not valid UTF-8");
    assert!(stderr.contains("Cannot access reference file"));
    assert!(stderr.contains("missing-reference.fna"));

    let _ = fs::remove_dir_all(temp_dir);
}

#[test]
fn direct_reference_paths_are_validated() {
    let exe = env!("CARGO_BIN_EXE_fasterANI");
    let query = fixture_path("Shigella_flexneri_2a_01.fna");

    let output = Command::new(exe)
        .args(["--reference", "missing-reference.fna", "--query", &query])
        .output()
        .expect("failed to launch fasterANI binary");

    assert!(
        !output.status.success(),
        "binary unexpectedly succeeded: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stderr = String::from_utf8(output.stderr).expect("stderr was not valid UTF-8");
    assert!(stderr.contains("Cannot access reference file"));
    assert!(stderr.contains("missing-reference.fna"));
}

#[test]
fn fasta_validation_rejects_non_file_path() {
    let exe = env!("CARGO_BIN_EXE_fasterANI");
    let temp_dir = temp_test_dir("non-file-reference");
    let query = fixture_path("Shigella_flexneri_2a_01.fna");

    let output = Command::new(exe)
        .args([
            "--reference",
            temp_dir.to_str().expect("utf-8 temp dir"),
            "--query",
            &query,
        ])
        .output()
        .expect("failed to launch fasterANI binary");

    assert!(
        !output.status.success(),
        "binary unexpectedly succeeded: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stderr = String::from_utf8(output.stderr).expect("stderr was not valid UTF-8");
    assert!(stderr.contains("Reference path is not a file"));

    let _ = fs::remove_dir_all(temp_dir);
}

#[test]
fn fasta_validation_rejects_files_of_100_bytes_or_less() {
    let exe = env!("CARGO_BIN_EXE_fasterANI");
    let reference = fixture_path("empty.fasta");
    let query = fixture_path("Shigella_flexneri_2a_01.fna");

    let output = Command::new(exe)
        .args(["--reference", &reference, "--query", &query])
        .output()
        .expect("failed to launch fasterANI binary");

    assert!(
        !output.status.success(),
        "binary unexpectedly succeeded: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stderr = String::from_utf8(output.stderr).expect("stderr was not valid UTF-8");
    assert!(stderr.contains("Reference file is too small"));
    assert!(stderr.contains("must be > 100"));
}

#[test]
fn skip_validation_bypasses_missing_direct_fasta_check() {
    let exe = env!("CARGO_BIN_EXE_fasterANI");
    let query = fixture_path("Shigella_flexneri_2a_01.fna");

    let output = Command::new(exe)
        .args([
            "--reference",
            "missing-reference.fna",
            "--query",
            &query,
            "--skip-validation",
        ])
        .output()
        .expect("failed to launch fasterANI binary");

    assert!(
        !output.status.success(),
        "binary unexpectedly succeeded with missing input"
    );
    let stderr = String::from_utf8(output.stderr).expect("stderr was not valid UTF-8");
    assert!(!stderr.contains("Cannot access reference file"));
    assert!(stderr.contains("# skip_validation = true  (CLI-only metadata)"));
}

#[test]
fn skip_validation_bypasses_missing_list_entry_check() {
    let exe = env!("CARGO_BIN_EXE_fasterANI");
    let temp_dir = temp_test_dir("skip-validation-list-entry");
    let reference_list = temp_dir.join("references.txt");
    let query = fixture_path("Shigella_flexneri_2a_01.fna");
    fs::write(&reference_list, "missing-reference.fna\n").expect("write reference list");

    let output = Command::new(exe)
        .args([
            "--reference-list",
            reference_list.to_str().expect("utf-8 reference list"),
            "--query",
            &query,
            "--skip-validation",
        ])
        .output()
        .expect("failed to launch fasterANI binary");

    assert!(
        !output.status.success(),
        "binary unexpectedly succeeded with missing list entry"
    );
    let stderr = String::from_utf8(output.stderr).expect("stderr was not valid UTF-8");
    assert!(!stderr.contains("Cannot access path 'missing-reference.fna' from list"));

    let _ = fs::remove_dir_all(temp_dir);
}

#[test]
fn release_documentation_covers_workflows_bounds_and_reproducibility() {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let readme = fs::read_to_string(manifest_dir.join("README.md")).expect("read README");
    let input_docs = fs::read_to_string(manifest_dir.join("docs/input-instructions.md"))
        .expect("read input instructions");

    for required in [
        "cargo install --path .",
        "Compare FASTA files directly",
        "Build an on-disk reference sketch",
        "Query that sketch later",
        "--minimizer-hash-seed",
    ] {
        assert!(readme.contains(required), "README omitted {required:?}");
    }
    for required in [
        "Complete Params File Example",
        "0..=4_294_967_295",
        "Valid range: 0..=100",
        "mphf_gamma = 10.0",
        "complete\neffective runtime configuration",
    ] {
        assert!(
            input_docs.contains(required),
            "input documentation omitted {required:?}"
        );
    }
}
