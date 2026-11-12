use std::fs;
use std::path::PathBuf;
use std::process::{Command, Stdio};

fn temp_dir(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("meta2bmap_cli_{}_{}", name, std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    dir
}

fn encode_base4(mut value: usize, length: usize) -> String {
    let mut encoded = vec![b'A'; length];
    for base in encoded.iter_mut().rev() {
        *base = b"ACGT"[value % 4];
        value /= 4;
    }
    String::from_utf8(encoded).unwrap()
}

fn bcgi_site(index: usize) -> String {
    format!(
        "{}CGA{}TGC{}",
        encode_base4(index, 10),
        encode_base4(index, 6),
        encode_base4(index + 1, 10)
    )
}

#[test]
fn extract_and_profile_minimal_native_2brad_sample() {
    let dir = temp_dir("end_to_end");
    let genome_path = dir.join("genome.fa");
    let reads_path = dir.join("sample.fq");
    let db_path = dir.join("run.syldb");
    let sample_path = dir.join("run.sylsp");
    let abundance_path = dir.join("abundance.tsv");

    let mut genome = String::from(">genome.fa\n");
    for index in 0..16 {
        genome.push_str(&bcgi_site(index));
        genome.push_str("AAAAAAAAAA");
    }
    fs::write(&genome_path, genome).unwrap();

    let mut reads = String::new();
    for index in 0..4 {
        let sequence = bcgi_site(index);
        reads.push_str(&format!(
            "@read{index}\n{sequence}\n+\n{}\n",
            "I".repeat(sequence.len())
        ));
    }
    fs::write(&reads_path, reads).unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_meta2bmap"))
        .args([
            "extract",
            "--threads",
            "1",
            "--genomes",
            genome_path.to_str().unwrap(),
            "--reads",
            reads_path.to_str().unwrap(),
            "--output",
            dir.to_str().unwrap(),
            "--sample-output-dir",
            dir.to_str().unwrap(),
            "--out-name",
            "run",
            "--enzyme",
            "BcgI",
        ])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "extract failed: {:?}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(db_path.is_file());
    assert!(sample_path.is_file());

    let profile_output = fs::File::create(dir.join("profile.tsv")).unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_meta2bmap"))
        .args([
            "profile",
            "--db-file",
            db_path.to_str().unwrap(),
            "--sample-file",
            sample_path.to_str().unwrap(),
            "--minimum-ani",
            "80",
            "--threads",
            "1",
            "--min-shared-tags",
            "1",
            "--min-tags-genome",
            "1",
            "--gscore-threshold",
            "0",
            "--tsv-name",
            abundance_path.to_str().unwrap(),
        ])
        .stdout(Stdio::from(profile_output))
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "profile failed: {:?}",
        output.stderr
    );
    let profile = fs::read_to_string(dir.join("profile.tsv")).unwrap();
    assert!(profile.contains("genome.fa"));
    assert!(abundance_path.is_file());

    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn rejects_invalid_sketch_parameters() {
    let dir = temp_dir("invalid_sketch");
    for args in [
        vec!["sketch", "--threads", "1", "--c-value", "0"],
        vec!["sketch", "--threads", "1", "--k-size", "32"],
        vec!["sketch", "--threads", "0"],
    ] {
        let output = Command::new(env!("CARGO_BIN_EXE_meta2bmap"))
            .args(&args)
            .output()
            .unwrap();
        assert!(
            !output.status.success(),
            "expected failure for {args:?}; stderr={:?}",
            String::from_utf8_lossy(&output.stderr)
        );
    }
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn reports_version() {
    let output = Command::new(env!("CARGO_BIN_EXE_meta2bmap"))
        .arg("--version")
        .output()
        .unwrap();
    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("meta2bmap"));
}
