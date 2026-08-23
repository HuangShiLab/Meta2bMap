// src/inspect.rs

use crate::cmdline::InspectArgs;
use anyhow::{Context, Result};
use bincode;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::Path;
// use regex::Regex;
use crate::constants::Hash;
use std::collections::HashMap;

#[derive(Serialize, Deserialize, Debug)]
struct InspectResult {
    file_type: String,
    file_name: String,
    enzyme: String,
    num_records: usize,
    total_tags: usize,
    unique_tags: usize,
    tag_frequency_stats: Vec<(Hash, usize)>,
    mean_read_length: Option<f64>,
    first_contig_name: Option<String>,
    genome_sources: Option<Vec<String>>,
    sample_sources: Option<Vec<SampleStats>>,
    // 每个样本对应的 tag -> 计数
    per_sample_tag_counts: Option<std::collections::HashMap<String, std::collections::HashMap<Hash, usize>>>,
    tag_lengths: Vec<usize>,
    tag_length_distribution: Vec<(usize, usize, f64)>,
    patterns: Vec<String>,
    genome_stats: Option<Vec<GenomeStats>>,
}

#[derive(Serialize, Deserialize, Debug)]
struct GenomeStats {
    source: String,
    num_records: usize,
    total_tags: usize,
    unique_tags: usize,
    tag_length_distribution: Vec<(usize, usize, f64)>,
}

#[derive(Serialize, Deserialize, Debug)]
struct SampleStats {
    source: String,
    num_records: usize,
    total_tags: usize,
    tag_length_distribution: Vec<(usize, usize, f64)>,
}

#[derive(Debug)]
struct TagMatrix {
    samples: Vec<String>,
    tags: Vec<Hash>,
    matrix: HashMap<(String, Hash), usize>,
}

pub fn inspect(args: InspectArgs) -> Result<()> {
    let mut writer = match args.out_file_name {
        Some(path) => Box::new(BufWriter::new(File::create(path)?)) as Box<dyn Write>,
        None => Box::new(BufWriter::new(std::io::stdout())) as Box<dyn Write>,
    };

    // 收集所有样本的标签信息用于生成矩阵
    let mut tag_matrix = TagMatrix {
        samples: Vec::new(),
        tags: Vec::new(),
        matrix: HashMap::new(),
    };

    for file in &args.files {
        match inspect_file(file) {
            Ok(result) => {
                // 输出文件信息
                writeln!(writer, "File Information:")?;
                writeln!(writer, "----------------")?;
                writeln!(writer, "File: {}", result.file_name)?;
                writeln!(writer, "Type: {}", result.file_type)?;
                // 样本有无与数量（若存在多个sample则按sample显示数量）
                match &result.sample_sources {
                    Some(samples) if !samples.is_empty() => {
                        writeln!(writer, "Samples: yes ({} samples)", samples.len())?;
                        writeln!(writer, "Per-sample total tags:")?;
                        for s in samples {
                            writeln!(writer, "  {}: {}", s.source, s.total_tags)?;
                        }
                    }
                    _ => {
                        writeln!(writer, "Samples: no")?;
                    }
                }
                writeln!(writer, "Total tags: {}", result.total_tags)?;
                writeln!(writer, "Unique tags: {}", result.unique_tags)?;
                
                writeln!(writer, "\nTag Frequency Statistics:")?;
                writeln!(writer, "------------------------")?;
                writeln!(writer, "{:<20} {:<10}", "Tag Hash", "Count")?;
                writeln!(writer, "{:-<30}", "")?;
                
                // 显示前20个最常见的tag
                let display_count = std::cmp::min(20, result.tag_frequency_stats.len());
                for (tag_hash, count) in result.tag_frequency_stats.iter().take(display_count) {
                    writeln!(writer, "{:<20} {:<10}", format!("{:016x}", tag_hash), count)?;
                }
                if result.tag_frequency_stats.len() > display_count {
                    writeln!(writer, "... and {} more unique tags", result.tag_frequency_stats.len() - display_count)?;
                }
                
                writeln!(writer, "\nDetected Enzyme Information:")?;
                writeln!(writer, "-------------------------")?;
                writeln!(writer, "Enzyme: {}", result.enzyme)?;
                writeln!(writer, "\nRecognition patterns:")?;
                for (i, pattern) in result.patterns.iter().enumerate() {
                    writeln!(writer, "  {}. {}", i + 1, pattern)?;
                }
                
                writeln!(writer, "\nTag Length Distribution:")?;
                writeln!(writer, "----------------------")?;
                writeln!(writer, "{:<10} {:<10}", "Length", "Count")?;
                writeln!(writer, "{:-<20}", "")?;
                
                for (length, count, _) in &result.tag_length_distribution {
                    writeln!(writer, "{:<10} {:<10}", length, count)?;
                }

                if let Some(sample_stats) = &result.sample_sources {
                    writeln!(writer, "\nSample-specific statistics:")?;
                    writeln!(writer, "------------------------")?;
                    for sample in sample_stats {
                        writeln!(writer, "\nSample: {}", sample.source)?;
                        writeln!(writer, "  Records: {}", sample.num_records)?;
                        writeln!(writer, "  Total tags: {}", sample.total_tags)?;
                        writeln!(writer, "  Tag length distribution:")?;
                        for (length, count, percentage) in &sample.tag_length_distribution {
                            writeln!(writer, "    Length {}: {} tags ({:.2}%)", length, count, percentage)?;
                        }
                    }
                }

                if let Some(stats) = &result.genome_stats {
                    writeln!(writer, "\nGenome-specific statistics:")?;
                    writeln!(writer, "------------------------")?;
                    for genome in stats {
                        // 提取基因组名称（去掉路径）
                        let genome_name = Path::new(&genome.source)
                            .file_stem()
                            .and_then(|s| s.to_str())
                            .unwrap_or(&genome.source);
                        writeln!(writer, "\nGenome: {}", genome_name)?;
                        writeln!(writer, "  Records: {}", genome.num_records)?;
                        writeln!(writer, "  Total tags: {}", genome.total_tags)?;
                        writeln!(writer, "  Unique tags: {}", genome.unique_tags)?;
                        writeln!(writer, "  Tag length distribution:")?;
                        for (length, count, _) in &genome.tag_length_distribution {
                            writeln!(writer, "    Length {}: {} tags", length, count)?;
                        }
                    }
                }
                writeln!(writer, "\n")?;

                // 收集标签矩阵信息
                collect_tag_matrix_data(&result, &mut tag_matrix);
            }
            Err(e) => {
                eprintln!("Failed to inspect {}: {}", file, e);
            }
        }
    }

    // 如果指定了输出路径，生成TSV矩阵
    if let Some(log_path) = &args.log_path {
        generate_tsv_matrix(&tag_matrix, log_path, &args.tsv_name)?;
    }

    Ok(())
}

fn inspect_file(file_path: &str) -> Result<InspectResult> {
    let path = Path::new(file_path);
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    match path.extension().and_then(|s| s.to_str()) {
        Some("db") | Some("syldb") => inspect_db(reader, file_path),
        Some("sp") | Some("sylsp") => inspect_sp(reader, file_path),
        _ => Err(anyhow::anyhow!("Unknown file extension, expected .db/.syldb or .sp/.sylsp")),
    }
}

fn inspect_db(reader: BufReader<File>, file_path: &str) -> Result<InspectResult> {
    let entries: Vec<crate::extract::SyldbEntry> = crate::extract::read_syldb(reader)
        .context("Failed to deserialize .db file")?;

    let mut tag_lengths = Vec::new();
    let mut genome_sources = std::collections::HashSet::new();
    let mut genome_stats = std::collections::HashMap::new();
    let mut tag_frequency = std::collections::HashMap::new();
    
    // 检查是否有unique标记信息
    let has_unique_marks = entries.iter().any(|entry| entry.tag_uniqueness.is_some());
    let mut total_unique_tags_marked = 0;

    for entry in &entries {
        let source = &entry.genome_source;
        genome_sources.insert(source.clone());
        
        let stats = genome_stats.entry(source.clone()).or_insert(GenomeStats {
            source: source.clone(),
            num_records: 0,
            total_tags: 0,
            unique_tags: 0,
            tag_length_distribution: Vec::new(),
        });
        
        stats.num_records += 1;
        stats.total_tags += entry.tags.len();
        
        // 处理tags和uniqueness信息
        for (i, tag) in entry.tags.iter().enumerate() {
            tag_lengths.push(8); // Hash is always 8 bytes (u64)
            *tag_frequency.entry(*tag).or_insert(0) += 1;
            
            // 如果有unique标记，统计unique tags
            if let Some(tag_uniqueness) = &entry.tag_uniqueness {
                if i < tag_uniqueness.len() && tag_uniqueness[i] {
                    stats.unique_tags += 1;
                    total_unique_tags_marked += 1;
                }
            }
        }
    }

    // 为每个基因组计算tag长度分布
    for stats in genome_stats.values_mut() {
        let mut lengths = Vec::new();
        for entry in &entries {
            if entry.genome_source == stats.source {
                for _tag in &entry.tags {
                    lengths.push(8); // Hash is always 8 bytes (u64)
                }
            }
        }
        stats.tag_length_distribution = calculate_tag_distribution(&lengths);
    }

    let distribution = calculate_tag_distribution(&tag_lengths);

    // 自描述数据库：从条目里读取构建时使用的酶；取众数（理论上所有条目相同）。
    let enzyme = entries.iter()
        .map(|e| e.enzyme.clone())
        .filter(|e| !e.is_empty())
        .max_by_key(|e| entries.iter().filter(|x| x.enzyme == *e).count())
        .unwrap_or_else(|| "unknown".to_string());
    let patterns = Vec::new();
    let _matched_count = 0usize;
    let _matched_ratio = 0.0;

    // 计算tag统计信息
    // 如果有unique标记，使用marked unique count；否则使用distinct tag count
    let unique_tags = if has_unique_marks {
        total_unique_tags_marked
    } else {
        tag_frequency.len()
    };
    
    let mut tag_frequency_stats: Vec<(Hash, usize)> = tag_frequency.into_iter().collect();
    tag_frequency_stats.sort_by(|a, b| b.1.cmp(&a.1)); // 按频率降序排序

    println!("Total tags: {}", tag_lengths.len());
    if has_unique_marks {
        println!("Unique tags: {} (marked as taxa-specific)", unique_tags);
    } else {
        println!("Unique tags: {} (distinct tag hashes)", unique_tags);
    }
    println!("Best matched enzyme: {}", enzyme);
    println!("Recognition patterns:");
    for (i, pattern) in patterns.iter().enumerate() {
        println!("  {}. {}", i + 1, pattern);
    }
    println!("Matched tags: {} ({:.2}%)", _matched_count, _matched_ratio * 100.0);

    Ok(InspectResult {
        file_type: "GenomeDatabase".to_string(),
        file_name: file_path.to_string(),
        enzyme,
        num_records: entries.len(),
        total_tags: tag_lengths.len(),
        unique_tags,
        tag_frequency_stats,
        mean_read_length: None,
        first_contig_name: entries.first().map(|e| e.sequence_id.clone()),
        genome_sources: if genome_sources.is_empty() {
            None
        } else {
            Some(genome_sources.into_iter().collect())
        },
        sample_sources: None,
        per_sample_tag_counts: None,
        tag_lengths,
        tag_length_distribution: distribution,
        patterns,
        genome_stats: Some(genome_stats.into_values().collect()),
    })
}

fn inspect_sp(reader: BufReader<File>, file_path: &str) -> Result<InspectResult> {
    let entries: Vec<crate::extract::SylspEntry> = bincode::deserialize_from(reader)
        .context("Failed to deserialize .sp file")?;

    let mut tag_lengths = Vec::new();
    let mut sample_stats = std::collections::HashMap::new();
    let mut tag_frequency = std::collections::HashMap::new();
    let mut per_sample_tag_counts: std::collections::HashMap<String, std::collections::HashMap<Hash, usize>> = std::collections::HashMap::new();

    for entry in &entries {
        tag_lengths.push(8); // Hash is always 8 bytes (u64)
        *tag_frequency.entry(entry.tag).or_insert(0) += 1;

        // 累积每个样本的 tag 计数
        let sample_entry = per_sample_tag_counts
            .entry(entry.sample_source.clone())
            .or_insert_with(std::collections::HashMap::new);
        *sample_entry.entry(entry.tag).or_insert(0) += 1;

        let stats = sample_stats.entry(entry.sample_source.clone()).or_insert(SampleStats {
            source: entry.sample_source.clone(),
            num_records: 0,
            total_tags: 0,
            tag_length_distribution: Vec::new(),
        });
        stats.num_records += 1;
        stats.total_tags += 1;
    }

    for stats in sample_stats.values_mut() {
        let mut sample_lengths = Vec::new();
        for entry in &entries {
            if entry.sample_source == stats.source {
                sample_lengths.push(8); // Hash is always 8 bytes (u64)
            }
        }
        stats.tag_length_distribution = calculate_tag_distribution(&sample_lengths);
    }

    let distribution = calculate_tag_distribution(&tag_lengths);
    // 样本文件自描述酶信息（新版格式存储；旧格式无此字段，显示 unknown）
    let enzyme = {
        let mut counts: HashMap<&str, usize> = HashMap::new();
        for e in &entries {
            if !e.enzyme.is_empty() {
                *counts.entry(e.enzyme.as_str()).or_insert(0) += 1;
            }
        }
        counts.into_iter().max_by_key(|(_, c)| *c).map(|(e, _)| e.to_string())
            .unwrap_or_else(|| "unknown".to_string())
    };
    let patterns = Vec::new();
    let _matched_count = 0usize;
    let _matched_ratio = 0.0;

    // 计算tag统计信息
    let unique_tags = tag_frequency.len();
    let mut tag_frequency_stats: Vec<(Hash, usize)> = tag_frequency.into_iter().collect();
    tag_frequency_stats.sort_by(|a, b| b.1.cmp(&a.1)); // 按频率降序排序

    Ok(InspectResult {
        file_type: "SampleProfile".to_string(),
        file_name: file_path.to_string(),
        enzyme,
        num_records: entries.len(),
        total_tags: tag_lengths.len(),
        unique_tags,
        tag_frequency_stats,
        mean_read_length: None,
        first_contig_name: entries.first().map(|e| e.sequence_id.clone()),
        genome_sources: None,
        sample_sources: Some(sample_stats.into_values().collect()),
        per_sample_tag_counts: Some(per_sample_tag_counts),
        tag_lengths,
        tag_length_distribution: distribution,
        patterns,
        genome_stats: None,
    })
}

// fn determine_enzyme_by_regex(tags: &[String]) -> (String, Vec<String>, usize, f64) {
//     if tags.is_empty() {
//         return ("unknown".to_string(), Vec::new(), 0, 0.0);
//     }
//     let mut best_enzyme = "unknown".to_string();
//     let mut best_patterns = Vec::new();
//     let mut best_count = 0;
//     let mut best_ratio = 0.0;
//     for (name, patterns) in crate::extract::ENZYME_DEFINITIONS {
//         let regexes: Vec<Regex> = patterns.iter().filter_map(|p| Regex::new(p).ok()).collect();
//         let mut count = 0;
//         for tag in tags {
//             if regexes.iter().any(|re| re.is_match(tag)) {
//                 count += 1;
//             }
//         }
//         let ratio = count as f64 / tags.len() as f64;
//         if count > best_count || (count == best_count && ratio > best_ratio) {
//             best_enzyme = name.to_string();
//             best_patterns = patterns.iter().map(|p| p.to_string()).collect();
//             best_count = count;
//             best_ratio = ratio;
//         }
//     }
//     (best_enzyme, best_patterns, best_count, best_ratio)
// }

fn calculate_tag_distribution(tag_lengths: &[usize]) -> Vec<(usize, usize, f64)> {
    use std::collections::HashMap;
    
    let mut length_counts: HashMap<usize, usize> = HashMap::new();
    for &length in tag_lengths {
        *length_counts.entry(length).or_insert(0) += 1;
    }

    let total = tag_lengths.len() as f64;
    let mut distribution: Vec<_> = length_counts
        .into_iter()
        .map(|(length, count)| {
            let percentage = (count as f64 / total) * 100.0;
            (length, count, percentage)
        })
        .collect();

    distribution.sort_by_key(|&(length, _, _)| length);
    distribution
}

fn collect_tag_matrix_data(result: &InspectResult, tag_matrix: &mut TagMatrix) {
    match result.file_type.as_str() {
        "SampleProfile" => {
            if let Some(per_sample) = &result.per_sample_tag_counts {
                for (sample_name, tag_counts) in per_sample {
                    if !tag_matrix.samples.contains(sample_name) {
                        tag_matrix.samples.push(sample_name.clone());
                    }
                    for (tag_hash, count) in tag_counts {
                        if !tag_matrix.tags.contains(tag_hash) {
                            tag_matrix.tags.push(*tag_hash);
                        }
                        let key = (sample_name.clone(), *tag_hash);
                        let entry = tag_matrix.matrix.entry(key).or_insert(0);
                        *entry += *count;
                    }
                }
            }
        }
        "GenomeDatabase" => {
            // 对于基因组数据库文件，使用文件名作为样本名
            let sample_name = Path::new(&result.file_name)
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("unknown")
                .to_string();
            
            if !tag_matrix.samples.contains(&sample_name) {
                tag_matrix.samples.push(sample_name.clone());
            }
            
            // 收集该数据库中的所有标签
            for (tag_hash, count) in &result.tag_frequency_stats {
                if !tag_matrix.tags.contains(tag_hash) {
                    tag_matrix.tags.push(*tag_hash);
                }
                let key = (sample_name.clone(), *tag_hash);
                let entry = tag_matrix.matrix.entry(key).or_insert(0);
                *entry += *count;
            }
        }
        _ => {}
    }
}

fn generate_tsv_matrix(tag_matrix: &TagMatrix, log_path: &str, tsv_name: &str) -> Result<()> {
    // 确保输出目录存在
    std::fs::create_dir_all(log_path)?;
    
    // 构建完整的文件路径
    let tsv_path = Path::new(log_path).join(tsv_name);
    let mut tsv_writer = BufWriter::new(File::create(&tsv_path)?);
    
    // 排序样本和标签以确保输出的一致性
    let mut sorted_samples = tag_matrix.samples.clone();
    sorted_samples.sort();
    let mut sorted_tags = tag_matrix.tags.clone();
    sorted_tags.sort();
    
    // 写入表头
    write!(tsv_writer, "Tag")?; // 目前保存的是哈希，无法还原原始序列，这里用16进制hash展示
    for sample in &sorted_samples {
        write!(tsv_writer, "\t{}", sample)?;
    }
    writeln!(tsv_writer)?;
    
    // 写入数据行
    for tag in &sorted_tags {
        write!(tsv_writer, "{:016x}", tag)?;
        for sample in &sorted_samples {
            let count = tag_matrix.matrix.get(&(sample.clone(), *tag)).unwrap_or(&0);
            write!(tsv_writer, "\t{}", count)?;
        }
        writeln!(tsv_writer)?;
    }
    
    println!("Tag count matrix saved to: {}", tsv_path.display());
    Ok(())
}
