use anyhow::{Context, Result};
use fxhash::FxHashMap;
use rayon::prelude::*;
use serde::{Serialize, Deserialize};
use std::{
    fs::File,
    io::{BufReader, BufWriter},
    path::Path,
};

use crate::cmdline::MarkArgs;
use crate::extract::SyldbEntry;
use crate::constants::Hash;

/// 包含unique标记统计信息的结构体
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct GenomeStats {
    pub genome_source: String,
    pub total_tags: usize,
    pub unique_tags: usize,
}

/// 标记unique tags的主函数
pub fn mark(args: MarkArgs) -> Result<()> {
    println!("开始标记unique tags...");
    
    // 读取.db文件
    let input_path = Path::new(&args.input_file);
    let db_entries = read_db_file(input_path)?;
    
    println!("已读取 {} 个db条目", db_entries.len());
    
    // 分析并标记unique tags
    let marked_entries = mark_unique_tags(db_entries)?;
    
    // 生成统计信息
    let stats = generate_statistics(&marked_entries);
    print_statistics(&stats);
    
    // 写回文件
    let output_path = if let Some(output) = args.output_file {
        Path::new(&output).to_path_buf()
    } else {
        input_path.to_path_buf()
    };
    
    write_db_file(&output_path, &marked_entries)?;
    
    println!("标记完成，已写入文件: {}", output_path.display());
    
    Ok(())
}

/// 读取.db文件
fn read_db_file(path: &Path) -> Result<Vec<SyldbEntry>> {
    let file = File::open(path)
        .context(format!("无法打开文件: {}", path.display()))?;
    let reader = BufReader::new(file);
    
    let entries: Vec<SyldbEntry> = crate::extract::read_syldb(reader)
        .context("无法反序列化db文件")?;
    
    Ok(entries)
}

/// 标记unique tags的核心逻辑
fn mark_unique_tags(mut entries: Vec<SyldbEntry>) -> Result<Vec<SyldbEntry>> {
    // genome_source 驻留（intern）为 u32 id，避免对每个 tag 克隆 String。
    // uniqueness 语义是「该 tag 只出现在 1 个基因组」，因此 tag 的状态只需记录
    // （首个基因组 id, 是否已见过第二个不同基因组），无需保存完整基因组集合。
    let mut genome_ids: FxHashMap<&str, u32> = FxHashMap::default();
    let mut next_gid = 0u32;
    // tag -> (first_genome_id, seen_multiple_genomes)
    let mut tag_state: FxHashMap<Hash, (u32, bool)> = FxHashMap::default();

    // 第一次遍历：收集所有tag和它们出现的基因组
    for entry in &entries {
        let gid = *genome_ids
            .entry(entry.genome_source.as_str())
            .or_insert_with(|| {
                let id = next_gid;
                next_gid += 1;
                id
            });
        for tag in &entry.tags {
            tag_state
                .entry(*tag)
                .and_modify(|(first, multi)| {
                    if *first != gid {
                        *multi = true;
                    }
                })
                .or_insert((gid, false));
        }
    }

    println!("总共找到 {} 个唯一tags", tag_state.len());

    // 计算unique tags数量
    let unique_tag_count = tag_state
        .values()
        .filter(|(_, multi)| !multi)
        .count();

    println!("其中 {} 个tags是unique的 ({:.2}%)",
        unique_tag_count,
        unique_tag_count as f64 / tag_state.len() as f64 * 100.0);

    // 第二次遍历：为每个entry标记其tags的uniqueness（tag_state 只读，可安全并行）
    entries.par_iter_mut().for_each(|entry| {
        let mut tag_uniqueness = Vec::with_capacity(entry.tags.len());

        for tag in &entry.tags {
            let is_unique = tag_state
                .get(tag)
                .map(|(_, multi)| !multi)
                .unwrap_or(false);
            tag_uniqueness.push(is_unique);
        }

        entry.tag_uniqueness = Some(tag_uniqueness);
    });

    Ok(entries)
}

/// 生成统计信息
fn generate_statistics(entries: &[SyldbEntry]) -> Vec<GenomeStats> {
    let mut genome_stats: FxHashMap<String, GenomeStats> = FxHashMap::default();
    
    for entry in entries {
        if let Some(tag_uniqueness) = &entry.tag_uniqueness {
            let stats = genome_stats
                .entry(entry.genome_source.clone())
                .or_insert_with(|| GenomeStats {
                    genome_source: entry.genome_source.clone(),
                    total_tags: 0,
                    unique_tags: 0,
                });
            
            stats.total_tags += entry.tags.len();
            stats.unique_tags += tag_uniqueness.iter().filter(|&&is_unique| is_unique).count();
        }
    }
    
    let mut stats: Vec<GenomeStats> = genome_stats.into_values().collect();
    stats.sort_by(|a, b| b.total_tags.cmp(&a.total_tags)); // 按total_tags降序排列
    
    stats
}

/// 打印统计信息
fn print_statistics(stats: &[GenomeStats]) {
    let total_tags: usize = stats.iter().map(|s| s.total_tags).sum();
    let total_unique_tags: usize = stats.iter().map(|s| s.unique_tags).sum();
    
    println!("\n=== Unique Tags 标记统计 ===");
    println!("Total tags: {}", total_tags);
    println!("Unique tags: {} ({:.2}%)", 
        total_unique_tags, 
        if total_tags > 0 { total_unique_tags as f64 / total_tags as f64 * 100.0 } else { 0.0 });
    
    println!("\nGenome-specific statistics:");
    println!("------------------------");
    
    for stat in stats {
        let unique_percentage = if stat.total_tags > 0 {
            stat.unique_tags as f64 / stat.total_tags as f64 * 100.0
        } else {
            0.0
        };
        
        // 提取基因组名称（去掉路径）
        let genome_name = Path::new(&stat.genome_source)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or(&stat.genome_source);
        
        println!("Genome: {}", genome_name);
        println!("  Total tags: {}", stat.total_tags);
        println!("  Unique tags: {} ({:.2}%)", stat.unique_tags, unique_percentage);
    }
}

/// 写入.db文件
fn write_db_file(path: &Path, entries: &[SyldbEntry]) -> Result<()> {
    let file = File::create(path)
        .context(format!("无法创建文件: {}", path.display()))?;
    let writer = BufWriter::new(file);
    
    crate::extract::write_syldb(writer, entries)
        .context("无法序列化db数据")?;
    
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::hash_bytes;
    
    #[test]
    fn test_mark_unique_tags() {
        let entries = vec![
            SyldbEntry {
                sequence_id: "seq1".to_string(),
                tags: vec![hash_bytes(b"ATGC"), hash_bytes(b"CGTA")],
                tag_lengths: vec![4, 4],
                genome_source: "genome_a.fa".to_string(),
                tag_uniqueness: None,
                tag_seqs: None,
                enzyme: "BcgI".to_string(),
                tag_positions: None,
                seq_len: 0,
            },
            SyldbEntry {
                sequence_id: "seq2".to_string(),
                tags: vec![hash_bytes(b"ATGC"), hash_bytes(b"TTTT")],
                tag_lengths: vec![4, 4],
                genome_source: "genome_b.fa".to_string(),
                tag_uniqueness: None,
                tag_seqs: None,
                enzyme: "BcgI".to_string(),
                tag_positions: None,
                seq_len: 0,
            },
        ];
        
        let marked_entries = mark_unique_tags(entries).unwrap();
        
        // ATGC在两个基因组中都出现，应该不是unique
        // CGTA只在genome_a中出现，应该是unique
        // TTTT只在genome_b中出现，应该是unique
        assert_eq!(marked_entries[0].tag_uniqueness.as_ref().unwrap()[0], false); // ATGC
        assert_eq!(marked_entries[0].tag_uniqueness.as_ref().unwrap()[1], true);  // CGTA
        assert_eq!(marked_entries[1].tag_uniqueness.as_ref().unwrap()[0], false); // ATGC
        assert_eq!(marked_entries[1].tag_uniqueness.as_ref().unwrap()[1], true);  // TTTT
    }
}
