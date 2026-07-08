//! End-to-end test that runs the compiled `fasterANI` binary against the bundled
//! test genomes and checks the emitted ANI line. This guards the public CLI
//! contract (arguments in, TSV out) the same way the README example does.

use std::{
    fs,
    io::Write,
    path::PathBuf,
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

fn assert_expected_test_data_result(stdout: &str) {
    let fields: Vec<&str> = stdout.trim_end().split('\t').collect();
    assert_eq!(fields.len(), 11, "unexpected result fields: {fields:?}");
    assert_eq!(fields[0], "assets/test-data/Shigella_flexneri_2a_01.fna");
    assert_eq!(
        fields[1],
        "assets/test-data/Escherichia_coli_str_K12_MG1655.fna"
    );
    assert_eq!(&fields[2..9], &["97.636", "0.807", "1608.00", "98.318", "2.500", "97.772", "97.500"]);

    for (index, field) in fields[9..].iter().enumerate() {
        field
            .parse::<f64>()
            .unwrap_or_else(|err| panic!("field {} was not numeric: {field:?}: {err}", index + 9));
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
    let sketch_prefix = sketch_prefix.to_str().expect("utf-8 sketch prefix");
    let mapping_stats_path = mapping_stats_path
        .to_str()
        .expect("utf-8 mapping stats path");
    let reference_one = "assets/test-data/Escherichia_coli_str_K12_MG1655.fna";
    let reference_two = "assets/test-data/Shigella_flexneri_2a_01.fna";

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
    assert!(!stderr.contains("FasterANI non-default runtime parameters"));
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
    assert!(!stderr.contains("FasterANI non-default runtime parameters"));
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
        .args([
            "--reference",
            &reference,
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
