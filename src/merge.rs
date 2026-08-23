use anyhow::{bail, Context, Result};
use fxhash::FxHashSet;
use std::{
    fs::File,
    io::{BufReader, BufWriter},
    path::Path,
};

use crate::cmdline::MergeArgs;
use crate::extract::SyldbEntry;

/// Merge multiple .syldb database files into one .syldb.
///
/// The current .syldb format (v2) is an 8-byte magic header followed by a
/// bincode-serialized `Vec<SyldbEntry>`; v1 files (no magic) are read via a
/// compatibility fallback in `read_syldb`. Merging is: deserialize each input,
/// concatenate, serialize once. Entries are per-fasta-record; batches are
/// assumed to contain disjoint genomes, so no deduplication is performed.
pub fn merge(args: MergeArgs) -> Result<()> {
    if args.inputs.len() < 2 {
        bail!("merge requires at least 2 input .syldb files");
    }

    let mut all_entries: Vec<SyldbEntry> = Vec::new();
    let mut enzymes: FxHashSet<String> = FxHashSet::default();
    let mut saw_empty_enzyme = false;
    let mut saw_tag_seqs = false;
    let mut saw_no_tag_seqs = false;

    for input in &args.inputs {
        let path = Path::new(input);
        let file = File::open(path)
            .with_context(|| format!("Failed to open input syldb file: {}", path.display()))?;
        let reader = BufReader::new(file);
        let entries: Vec<SyldbEntry> = crate::extract::read_syldb(reader)
            .with_context(|| format!("Failed to deserialize syldb file: {}", path.display()))?;

        for e in &entries {
            if e.enzyme.is_empty() {
                saw_empty_enzyme = true;
            } else {
                enzymes.insert(e.enzyme.clone());
            }
            if e.tag_seqs.is_some() {
                saw_tag_seqs = true;
            } else {
                saw_no_tag_seqs = true;
            }
        }

        eprintln!(
            "{}: {} entries",
            path.display(),
            entries.len()
        );
        all_entries.extend(entries);
    }

    if enzymes.len() > 1 {
        bail!(
            "Input databases were built with different enzymes: {:?}. Refusing to merge.",
            enzymes
        );
    }
    if saw_empty_enzyme && !enzymes.is_empty() {
        eprintln!(
            "Warning: some input databases lack enzyme metadata (built by an older version); assuming they match enzyme '{}'.",
            enzymes.iter().next().unwrap()
        );
    }
    if saw_tag_seqs && saw_no_tag_seqs {
        eprintln!(
            "Warning: mixing databases with and without tag sequences; --mismatch 1 will only work for entries that carry tag_seqs."
        );
    }

    let out_path = Path::new(&args.output);
    let out_file = File::create(out_path)
        .with_context(|| format!("Failed to create output syldb file: {}", out_path.display()))?;
    let writer = BufWriter::new(out_file);
    crate::extract::write_syldb(writer, &all_entries)
        .context("Failed to serialize merged syldb data")?;

    eprintln!(
        "Merged {} input files ({} total entries) -> {}",
        args.inputs.len(),
        all_entries.len(),
        out_path.display()
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::Hash;

    fn make_entry(id: &str, enzyme: &str, with_seqs: bool) -> SyldbEntry {
        let tags: Vec<Hash> = vec![1, 2, 3];
        SyldbEntry {
            sequence_id: id.to_string(),
            tags: tags.clone(),
            tag_lengths: vec![32, 32, 32],
            genome_source: format!("{}.fa", id),
            tag_uniqueness: None,
            tag_seqs: if with_seqs {
                Some(tags.iter().map(|t| t.to_le_bytes().to_vec()).collect())
            } else {
                None
            },
            enzyme: enzyme.to_string(),
            tag_positions: None,
            seq_len: 0,
        }
    }

    #[test]
    fn test_merge_roundtrip_preserves_entries() {
        let dir = std::env::temp_dir();
        let p1 = dir.join("m2b_merge_test_a.syldb");
        let p2 = dir.join("m2b_merge_test_b.syldb");
        let pout = dir.join("m2b_merge_test_out.syldb");

        let entries_a = vec![make_entry("g1", "BcgI", true), make_entry("g2", "BcgI", true)];
        let entries_b = vec![make_entry("g3", "BcgI", false)];

        crate::extract::write_syldb(BufWriter::new(File::create(&p1).unwrap()), &entries_a).unwrap();
        crate::extract::write_syldb(BufWriter::new(File::create(&p2).unwrap()), &entries_b).unwrap();

        let args = MergeArgs {
            inputs: vec![
                p1.to_string_lossy().to_string(),
                p2.to_string_lossy().to_string(),
            ],
            output: pout.to_string_lossy().to_string(),
        };
        merge(args).unwrap();

        let merged: Vec<SyldbEntry> =
            crate::extract::read_syldb(BufReader::new(File::open(&pout).unwrap())).unwrap();
        assert_eq!(merged.len(), 3);
        assert_eq!(merged[0].sequence_id, "g1");
        assert_eq!(merged[1].sequence_id, "g2");
        assert_eq!(merged[2].sequence_id, "g3");
        assert!(merged[0].tag_seqs.is_some());
        assert!(merged[2].tag_seqs.is_none());
        assert_eq!(merged[2].tags, vec![1, 2, 3]);

        let _ = std::fs::remove_file(&p1);
        let _ = std::fs::remove_file(&p2);
        let _ = std::fs::remove_file(&pout);
    }

    #[test]
    fn test_merge_rejects_mismatched_enzymes() {
        let dir = std::env::temp_dir();
        let p1 = dir.join("m2b_merge_test_e1.syldb");
        let p2 = dir.join("m2b_merge_test_e2.syldb");
        let pout = dir.join("m2b_merge_test_eout.syldb");

        crate::extract::write_syldb(
            BufWriter::new(File::create(&p1).unwrap()),
            &vec![make_entry("g1", "BcgI", true)],
        )
        .unwrap();
        crate::extract::write_syldb(
            BufWriter::new(File::create(&p2).unwrap()),
            &vec![make_entry("g2", "BslFI", true)],
        )
        .unwrap();

        let args = MergeArgs {
            inputs: vec![
                p1.to_string_lossy().to_string(),
                p2.to_string_lossy().to_string(),
            ],
            output: pout.to_string_lossy().to_string(),
        };
        assert!(merge(args).is_err());

        let _ = std::fs::remove_file(&p1);
        let _ = std::fs::remove_file(&p2);
        let _ = std::fs::remove_file(&pout);
    }
}
