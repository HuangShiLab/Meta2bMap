use crate::cmdline::QueryArgs;
use anyhow::{Result, Context};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write, BufReader};
use bincode;
use crate::constants::Hash;

pub fn query(args: QueryArgs) -> Result<()> {
    let mut writer = match args.out_file_name {
        Some(path) => Box::new(BufWriter::new(File::create(path)?)) as Box<dyn Write>,
        None => Box::new(BufWriter::new(std::io::stdout())) as Box<dyn Write>,
    };

    // 读取数据库文件
    let db_file = File::open(&args.db_file)?;
    let db_reader = BufReader::new(db_file);
    let db_entries: Vec<crate::extract::SyldbEntry> = crate::extract::read_syldb(db_reader)
        .context("Failed to deserialize database file")?;

    println!("Debug: Found {} entries in database file", db_entries.len());

    // 读取查询文件
    let query_file = File::open(&args.query_file)?;
    let query_reader = BufReader::new(query_file);
    let query_entries: Vec<crate::extract::SylspEntry> = bincode::deserialize_from(query_reader)
        .context("Failed to deserialize query file")?;

    println!("Debug: Found {} entries in query file", query_entries.len());

    // 创建基因组标签映射
    let mut genome_tags: HashMap<String, HashSet<Hash>> = HashMap::new();
    for entry in &db_entries {
        let tags = genome_tags.entry(entry.genome_source.clone())
            .or_insert_with(HashSet::new);
        tags.extend(entry.tags.iter().cloned());
    }

    println!("Debug: Found {} unique genomes in database", genome_tags.len());
    for (genome, tags) in &genome_tags {
        println!("Debug: Genome {} has {} tags", genome, tags.len());
    }

    // 收集查询标签
    let query_tags: HashSet<Hash> = query_entries.iter()
        .map(|e| e.tag.clone())
        .collect();

    println!("Debug: Found {} unique tags in query", query_tags.len());

    // 计算每个基因组的 ANI
    let mut results = Vec::new();
    
    for (genome, tags) in &genome_tags {
        let common_tags: HashSet<_> = query_tags.intersection(tags).collect();
        let ani = if !tags.is_empty() {
            (common_tags.len() as f64 / tags.len() as f64) * 100.0
        } else {
            0.0
        };

        println!("Debug: Genome {} - ANI: {:.2}%, Common tags: {}, Total tags: {}", 
            genome, ani, common_tags.len(), tags.len());

        // 只处理 ANI 大于最小阈值的基因组
        if ani >= args.minimum_ani {
            results.push((genome.clone(), ani, common_tags.len(), tags.len()));
        }
    }

    println!("Debug: Found {} genomes above minimum ANI threshold", results.len());

    // 按 ANI 降序排序
    results.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

    // 输出结果
    writeln!(writer, "Query Results:")?;
    writeln!(writer, "-------------")?;
    writeln!(writer, "Query file: {}", args.query_file)?;
    writeln!(writer, "Database file: {}", args.db_file)?;
    writeln!(writer, "\nGenome matches:")?;
    writeln!(writer, "{:<50} {:<10} {:<15} {:<15}", 
        "Genome", "ANI (%)", "Common Tags", "Total Tags")?;
    writeln!(writer, "{:-<90}", "")?;

    for (genome, ani, common_tags, total_tags) in &results {
        writeln!(writer, "{:<50} {:<10.2} {:<15} {:<15}", 
            genome, ani, common_tags, total_tags)?;
    }

    Ok(())
} 