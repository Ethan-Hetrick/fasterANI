//! End-to-end test that runs the compiled `fasterANI` binary against the bundled
//! test genomes and checks the emitted ANI line. This guards the public CLI
//! contract (arguments in, TSV out) the same way the README example does.

use std::process::Command;

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
assets/test-data/Escherichia_coli_str_K12_MG1655.fna\t97.636\t1297.00\t1608.00\n";
    assert_eq!(stdout, expected);
}
