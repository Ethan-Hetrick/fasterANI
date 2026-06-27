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
    let expected = "assets/test-data/Shigella_flexneri_2a_01.fna\t\
assets/test-data/Escherichia_coli_str_K12_MG1655.fna\t97.636\t1297.00\t1608.00\t98.318\t2.500\t97.772\t97.500\n";
    assert_eq!(stdout, expected);
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
    let expected = "assets/test-data/Shigella_flexneri_2a_01.fna\t\
assets/test-data/Escherichia_coli_str_K12_MG1655.fna\t97.636\t1297.00\t1608.00\t98.318\t2.500\t97.772\t97.500\n";
    assert_eq!(stdout, expected);
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
