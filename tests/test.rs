// tests/integration_test.rs
use assert_cmd::Command;
use std::fs;
use std::path::PathBuf;

// Global setup - runs once before any tests
#[ctor::ctor]
fn setup() {
    let output_dir = PathBuf::from("tests/output");
    // Clean up any existing output directory
    if output_dir.exists() {
        fs::remove_dir_all(&output_dir).ok();
    }
    fs::create_dir_all(&output_dir).expect("Failed to create output directory");
}

// Global cleanup - runs once after all tests complete
#[ctor::dtor]
fn cleanup() {
    let output_dir = PathBuf::from("tests/output");
    if output_dir.exists() {
        fs::remove_dir_all(&output_dir).ok();
    }
}

fn compare_files(output_path: &str, expected_path: &str) {
    let output = fs::read_to_string(output_path)
        .unwrap_or_else(|_| panic!("Failed to read output file: {}", output_path));
    let expected = fs::read_to_string(expected_path)
        .unwrap_or_else(|_| panic!("Failed to read expected file: {}", expected_path));
    
    assert_eq!(
        output, expected,
        "File mismatch:\nOutput: {}\nExpected: {}",
        output_path, expected_path
    );
}

fn check_stderr(stderr: &[u8]) {
    let stderr_str = String::from_utf8_lossy(stderr);
    let expected = fs::read_to_string("tests/expected/test_all.stderr")
        .expect("Failed to read tests/expected/test_all.stderr");
    
    assert_eq!(
        stderr_str, expected,
        "Stderr mismatch:\nActual stderr:\n{}\nExpected stderr:\n{}",
        stderr_str, expected
    );
}

#[test]
fn test_plain() {
    let suffix = "plain";
    let mut cmd = Command::cargo_bin("maple_to_vcf").unwrap();
    
    let assert = cmd.arg("--input-directory")
        .arg("tests/input")
        .arg("--ref-fasta")
        .arg("tests/input/ref.fasta")
        .arg("--output")
        .arg(format!("tests/output/test_{}.vcf", suffix))
        .assert()
        .success();
    
    check_stderr(&assert.get_output().stderr);
    
    compare_files(
        &format!("tests/output/test_{}.vcf", suffix),
        &format!("tests/expected/test_{}.vcf", suffix),
    );
}

#[test]
fn test_exclude() {
    let suffix = "exclude";
    let mut cmd = Command::cargo_bin("maple_to_vcf").unwrap();
    
    let assert = cmd.arg("--input-directory")
        .arg("tests/input")
        .arg("--ref-fasta")
        .arg("tests/input/ref.fasta")
        .arg("--exclude-file")
        .arg("tests/input/exclude.txt")
        .arg("--output")
        .arg(format!("tests/output/test_{}.vcf", suffix))
        .assert()
        .success();
    
    check_stderr(&assert.get_output().stderr);
    
    compare_files(
        &format!("tests/output/test_{}.vcf", suffix),
        &format!("tests/expected/test_{}.vcf", suffix),
    );
}

#[test]
fn test_min_maf() {
    let suffix = "min_maf";
    let mut cmd = Command::cargo_bin("maple_to_vcf").unwrap();
    
    let assert = cmd.arg("--input-directory")
        .arg("tests/input")
        .arg("--ref-fasta")
        .arg("tests/input/ref.fasta")
        .arg("--min-maf")
        .arg("0.1")
        .arg("--output")
        .arg(format!("tests/output/test_{}.vcf", suffix))
        .assert()
        .success();
    
    check_stderr(&assert.get_output().stderr);
    
    compare_files(
        &format!("tests/output/test_{}.vcf", suffix),
        &format!("tests/expected/test_{}.vcf", suffix),
    );
}

#[test]
fn test_min_maf_exclude() {
    let suffix = "min_maf_exclude";
    let mut cmd = Command::cargo_bin("maple_to_vcf").unwrap();
    
    let assert = cmd.arg("--input-directory")
        .arg("tests/input")
        .arg("--ref-fasta")
        .arg("tests/input/ref.fasta")
        .arg("--exclude-file")
        .arg("tests/input/exclude.txt")
        .arg("--min-maf")
        .arg("0.1")
        .arg("--output")
        .arg(format!("tests/output/test_{}.vcf", suffix))
        .assert()
        .success();
    
    check_stderr(&assert.get_output().stderr);
    
    compare_files(
        &format!("tests/output/test_{}.vcf", suffix),
        &format!("tests/expected/test_{}.vcf", suffix),
    );
}

#[test]
fn test_high_min_maf() {
    let suffix = "high_min_maf";
    let mut cmd = Command::cargo_bin("maple_to_vcf").unwrap();
    
    let assert = cmd.arg("--input-directory")
        .arg("tests/input")
        .arg("--ref-fasta")
        .arg("tests/input/ref.fasta")
        .arg("--min-maf")
        .arg("0.35")
        .arg("--output")
        .arg(format!("tests/output/test_{}.vcf", suffix))
        .assert()
        .success();
    
    check_stderr(&assert.get_output().stderr);
    
    compare_files(
        &format!("tests/output/test_{}.vcf", suffix),
        &format!("tests/expected/test_{}.vcf", suffix),
    );
}

#[test]
fn test_high_min_maf_exclude() {
    let suffix = "high_min_maf_exclude";
    let mut cmd = Command::cargo_bin("maple_to_vcf").unwrap();
    
    let assert = cmd.arg("--input-directory")
        .arg("tests/input")
        .arg("--ref-fasta")
        .arg("tests/input/ref.fasta")
        .arg("--min-maf")
        .arg("0.35")
        .arg("--exclude-file")
        .arg("tests/input/exclude.txt")
        .arg("--output")
        .arg(format!("tests/output/test_{}.vcf", suffix))
        .assert()
        .success();
    
    check_stderr(&assert.get_output().stderr);
    
    compare_files(
        &format!("tests/output/test_{}.vcf", suffix),
        &format!("tests/expected/test_{}.vcf", suffix),
    );
}

#[test]
fn test_min_non_n() {
    let suffix = "min_non_N";
    let mut cmd = Command::cargo_bin("maple_to_vcf").unwrap();
    
    let assert = cmd.arg("--input-directory")
        .arg("tests/input")
        .arg("--ref-fasta")
        .arg("tests/input/ref.fasta")
        .arg("--min-non-N")
        .arg("0.75")
        .arg("--output")
        .arg(format!("tests/output/test_{}.vcf", suffix))
        .assert()
        .success();
    
    check_stderr(&assert.get_output().stderr);
    
    compare_files(
        &format!("tests/output/test_{}.vcf", suffix),
        &format!("tests/expected/test_{}.vcf", suffix),
    );
}

#[test]
fn test_min_non_n_exclude() {
    let suffix = "min_non_N_exclude";
    let mut cmd = Command::cargo_bin("maple_to_vcf").unwrap();
    
    let assert = cmd.arg("--input-directory")
        .arg("tests/input")
        .arg("--ref-fasta")
        .arg("tests/input/ref.fasta")
        .arg("--exclude-file")
        .arg("tests/input/exclude.txt")
        .arg("--min-non-N")
        .arg("0.75")
        .arg("--output")
        .arg(format!("tests/output/test_{}.vcf", suffix))
        .assert()
        .success();
    
    check_stderr(&assert.get_output().stderr);
    
    compare_files(
        &format!("tests/output/test_{}.vcf", suffix),
        &format!("tests/expected/test_{}.vcf", suffix),
    );
}
