use crate::cmdline::{ContainArgs, ProfileArgs};
use crate::constants::{hash_bytes, Hash};
use anyhow::{anyhow, Context, Result};
use fxhash::FxHashMap;
use memory_stats::memory_stats;
use rayon::prelude::*;
use std::collections::HashMap;
use std::collections::HashSet;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Write};
use std::path::PathBuf;
use std::sync::Arc;
use std::sync::Mutex;
use std::thread;
use std::time::Duration;

use crate::extract::GenomeSketch;
pub use crate::extract::{
    ani_from_containment_one_mismatch, ani_from_containment_one_mismatch_err,
    one_mismatch_canonical_hashes, SyldbEntry, SylspEntry, TagHash,
};
use crate::sketch::SequencesSketch;

// 内存监控和限制功能 - 采用 sylph 的真正实现
// pub fn check_vram_and_block(max_ram: usize, file: &str) {
//     if let Some(usage) = memory_stats() {
//         let mut gb_usage_curr = usage.virtual_mem as f64 / 1_000_000_000 as f64;
//         if (max_ram as f64) < gb_usage_curr {
//             eprintln!("Max memory reached. Blocking processing for {}. Curr memory {:.2} GB, max mem {} GB",
//                      file, gb_usage_curr, max_ram);
//         }
//         while (max_ram as f64) < gb_usage_curr {
//             let five_second = Duration::from_secs(1);
//             thread::sleep(five_second);
//             if let Some(usage) = memory_stats() {
//                 gb_usage_curr = usage.virtual_mem as f64 / 1_000_000_000 as f64;
//                 if (max_ram as f64) >= gb_usage_curr {
//                     eprintln!("Processing for {} freed", file);
//                 }
//             } else {
//                 break;
//             }
//         }
//     }
// }

// 定义分类学信息结构体
#[derive(Debug, Clone, Default)]
pub struct TaxonomyInfo {
    pub kingdom: String,    // d__
    pub phylum: String,     // p__
    pub class: String,      // c__
    pub order: String,      // o__
    pub family: String,     // f__
    pub genus: String,      // g__
    pub species: String,    // s__
}

impl TaxonomyInfo {
    pub fn from_gtdb_string(gtdb_str: &str) -> Result<Self> {
        let mut taxonomy = TaxonomyInfo::default();
        
        for part in gtdb_str.split(';') {
            let part = part.trim();
            if part.starts_with("d__") {
                taxonomy.kingdom = part[3..].to_string();
            } else if part.starts_with("p__") {
                taxonomy.phylum = part[3..].to_string();
            } else if part.starts_with("c__") {
                taxonomy.class = part[3..].to_string();
            } else if part.starts_with("o__") {
                taxonomy.order = part[3..].to_string();
            } else if part.starts_with("f__") {
                taxonomy.family = part[3..].to_string();
            } else if part.starts_with("g__") {
                taxonomy.genus = part[3..].to_string();
            } else if part.starts_with("s__") {
                taxonomy.species = part[3..].to_string();
            }
        }
        
        Ok(taxonomy)
    }
    
    pub fn get_species_key(&self) -> String {
        format!("{}|{}|{}|{}|{}|{}|{}", 
                self.kingdom, self.phylum, self.class, 
                self.order, self.family, self.genus, self.species)
    }
}

// 物种级别的丰度结果 - 使用 Arc 共享 TaxonomyInfo
#[derive(Debug, Clone)]
pub struct SpeciesAbundanceResult {
    pub taxonomy: Arc<TaxonomyInfo>,
    pub sample_abundances: FxHashMap<String, f64>,  // sample_id -> abundance
    pub total_tags: usize,
    pub genome_count: usize,  // 该物种包含的genome数量
    pub reads_count: usize,   // 该物种的总 reads 数 (S)
    pub gscore: f64,          // G-score = sqrt(reads_count * total_tags)
}

// 定义比对结果结构
#[derive(Debug, Clone)]
pub struct QueryResult {
    pub sample_file: String,
    pub genome_file: String,
    pub adjusted_ani: f64,
    pub eff_cov: f64,
    pub ani_percentile: (f64, f64),
    pub eff_lambda: f64,
    pub lambda_percentile: (f64, f64),
    pub median_cov: f64,
    pub mean_cov_geq1: f64,
    pub containment_ind: String,
    pub naive_ani: f64,
    pub contig_name: String,
    pub ref_tags: usize,
    pub shared_tags: usize,
    pub query_tags: usize,
    pub taxonomic_abundance: f64,
    pub sequence_abundance: f64,
}

// 新增基因组级别的结果结构体
#[derive(Debug, Clone)]
pub struct GenomeProfileResult {
    pub genome_id: String,
    pub sample_id: String,  // 这个字段将存储实际的样本来源（如 sample1, sample4 等）
    pub file_path: String,  // 新增字段，存储文件路径
    pub adjusted_ani: f64,
    pub taxonomic_abundance: f64,
    pub sequence_abundance: f64,
    pub common_tags: usize,
    pub total_tags: usize,
    pub eff_cov: f64,
}

// 新增：k-mer重新分配相关的结构体和函数

// Winner table条目结构 - 对应sylph中的 (f64, &'a GenomeSketch, bool)
// genome 用 u32 索引（通过 GenomeInterner 驻留），而不是每个 tag 都存一份 String，
// winner map 可能有上千万条目，String→u32 大幅降低内存并把比较变成整数比较。
#[derive(Debug, Clone)]
struct WinnerTableEntry {
    pub ani: f64,
    pub genome_idx: u32,
    pub was_reassigned: bool,
}

/// 把 genome_id 字符串驻留成 u32 索引，供 winner table 复用。
struct GenomeInterner {
    ids: Vec<String>,
    idx: FxHashMap<String, u32>,
}

impl GenomeInterner {
    fn new() -> Self {
        Self { ids: Vec::new(), idx: FxHashMap::default() }
    }
    fn intern(&mut self, id: &str) -> u32 {
        if let Some(&i) = self.idx.get(id) {
            return i;
        }
        let i = self.ids.len() as u32;
        self.ids.push(id.to_string());
        self.idx.insert(id.to_string(), i);
        i
    }
    fn get(&self, id: &str) -> Option<u32> {
        self.idx.get(id).copied()
    }
    fn resolve(&self, i: u32) -> &str {
        &self.ids[i as usize]
    }
}

// 重新分配统计信息
#[derive(Debug, Clone)]
struct ReassignmentStats {
    pub tags_lost: usize,
    pub total_tags: usize,
    pub reassignment_ratio: f64,
}

// 构建winner table - 借鉴sylph的高效实现
fn build_winner_table<'a>(
    results: &'a [QueryResult],
    db_entries: &'a [SyldbEntry],
    log_reassign: bool
) -> (FxHashMap<Hash, WinnerTableEntry>, GenomeInterner) {
    eprintln!("Building winner table for {} results and {} database entries", results.len(), db_entries.len());

    let mut tag_to_genome_map: FxHashMap<Hash, WinnerTableEntry> = FxHashMap::default();
    let mut interner = GenomeInterner::new();

    // 关键优化1：构建sequence_id到db_entry的直接映射，避免O(G)查找
    let mut seq_id_to_entry: FxHashMap<String, &SyldbEntry> = FxHashMap::default();
    for entry in db_entries {
        seq_id_to_entry.insert(entry.sequence_id.clone(), entry);
    }

    eprintln!("Built seq_id_to_entry mapping with {} entries", seq_id_to_entry.len());

    // 关键优化2：借鉴sylph的直接遍历方式，避免复杂的嵌套查找
    for res in results.iter() {
        // O(1)查找，而不是O(G)的find操作
        if let Some(db_entry) = seq_id_to_entry.get(&res.contig_name) {
            // 每个基因组只驻留一次，tag 循环里只用 u32，避免逐 tag 的 String 分配
            let genome_idx = interner.intern(extract_genome_id_from_path(&db_entry.genome_source));

            // 借鉴sylph的简洁遍历方式
            for tag in &db_entry.tags {
                let entry = tag_to_genome_map.entry(*tag).or_insert_with(|| {
                    WinnerTableEntry {
                        ani: res.adjusted_ani,
                        genome_idx,
                        was_reassigned: false,
                    }
                });

                // 关键：选择ANI最高的基因组作为该tag的"赢家"
                if res.adjusted_ani > entry.ani {
                    *entry = WinnerTableEntry {
                        ani: res.adjusted_ani,
                        genome_idx,
                        was_reassigned: true,
                    };
                }
            }
        } else {
            eprintln!("Warning: No database entry found for contig {}", res.contig_name);
        }
    }

    eprintln!("Winner table built with {} unique tags", tag_to_genome_map.len());

    // 记录重新分配日志（借鉴sylph的简洁实现）
    if log_reassign {
        eprintln!("------------- Logging tag reassignments -----------------");
        let mut genome_to_index: FxHashMap<String, usize> = FxHashMap::default();
        for (i, res) in results.iter().enumerate() {
            eprintln!("Index\t{}\t{}\t{}", i, res.genome_file, res.contig_name);
            genome_to_index.insert(res.genome_file.clone(), i);
        }

        // 关键优化3：借鉴sylph的简洁并行计算，避免复杂的映射查找
        (0..results.len()).into_par_iter().for_each(|i| {
            let res = &results[i];
            let mut reassign_edge_map: FxHashMap<(usize, usize), usize> = FxHashMap::default();

            // 使用优化后的seq_id_to_entry映射
            if let Some(db_entry) = seq_id_to_entry.get(&res.contig_name) {
                for tag in &db_entry.tags {
                    if let Some(winner_entry) = tag_to_genome_map.get(tag) {
                        let winner_genome = interner.resolve(winner_entry.genome_idx);
                        if winner_genome != res.genome_file {
                            if let Some(&winner_index) = genome_to_index.get(winner_genome) {
                                let edge_count = reassign_edge_map.entry((winner_index, i)).or_insert(0);
                                *edge_count += 1;
                            }
                        }
                    }
                }
            }

            // 直接输出，避免收集到向量中
            for ((from_idx, to_idx), count) in reassign_edge_map {
                if count > 10 {
                    eprintln!("{}->{}\t{}\ttags reassigned", from_idx, to_idx, count);
                }
            }
        });
    }

    // 添加调试信息：统计重新分配的情况
    let mut reassigned_tags = 0;
    let mut total_tags = 0;
    for entry in tag_to_genome_map.values() {
        total_tags += 1;
        if entry.was_reassigned {
            reassigned_tags += 1;
        }
    }
    eprintln!("Reassignment statistics: {}/{} tags were reassigned ({:.2}%)",
              reassigned_tags, total_tags,
              if total_tags > 0 { reassigned_tags as f64 / total_tags as f64 * 100.0 } else { 0.0 });

    // tag 重叠分析纯属诊断信息，却要对所有基因组的所有 tag 再建一次 FxHashMap，
    // 代价高。仅在 debug 日志开启时执行。
    if log::log_enabled!(log::Level::Debug) {
        let mut tag_counts: FxHashMap<Hash, usize> = FxHashMap::default();
        for db_entry in db_entries {
            for tag in &db_entry.tags {
                *tag_counts.entry(*tag).or_insert(0) += 1;
            }
        }

        let overlapping_tags = tag_counts.values().filter(|&&count| count > 1).count();
        let total_unique_tags = tag_counts.len();
        eprintln!("Tag overlap analysis: {}/{} unique tags are shared between genomes ({:.2}%)",
                  overlapping_tags, total_unique_tags,
                  if total_unique_tags > 0 { overlapping_tags as f64 / total_unique_tags as f64 * 100.0 } else { 0.0 });

        if overlapping_tags > 0 {
            eprintln!("Shared tags found! This should enable reassignment.");
        } else {
            eprintln!("No shared tags found. Reassignment cannot work without overlapping tags.");
        }
    }

    (tag_to_genome_map, interner)
}

// Per-sample winner-table recalculation — the core of recalculate_with_winner_table.
// Shared by the 2bRAD path (counts built from SylspEntry lists) and the sketch path
// (counts borrowed directly from SequencesSketch::kmer_counts).
#[allow(clippy::too_many_arguments)]
fn recalculate_one_sample_with_winner_table(
    sample_source: &str,
    sample_counts: &FxHashMap<Hash, u32>,
    total_sample_tags: usize,
    db_entries: &[SyldbEntry],
    winner_map: &FxHashMap<Hash, WinnerTableEntry>,
    interner: &GenomeInterner,
    min_ani: f64,
    log_reassign: bool,
    mismatch: usize,
    min_shared_tags: usize,
    min_tags_genome: usize,
    read_error_rate: Option<f64>,
    no_error_correction: bool,
    initial_by_contig: &FxHashMap<&str, &QueryResult>,
    reassign_protection: bool,
) -> Vec<QueryResult> {
    let error_rate = sample_error_rate(
        db_entries,
        sample_counts,
        mismatch,
        read_error_rate,
        no_error_correction,
    );
    let mut results = Vec::new();

    // 为每个基因组（已聚合）计算重新分配后的结果
    for db_entry in db_entries {
        // 最小标签数过滤
        if db_entry.tags.len() < min_tags_genome {
            continue;
        }

        let genome_id = extract_genome_id_from_path(&db_entry.genome_source);
        // 当前基因组在 winner table 中的索引（若从未赢得任何 tag 则为 None）
        let my_idx = interner.get(genome_id);
        let mut covs: Vec<u32> = Vec::new();
        let mut tags_lost_count = 0;

        // 只统计被 winner table 判给当前基因组的 tag 的覆盖度
        for (i, tag) in db_entry.tags.iter().enumerate() {
            let seq = db_entry.tag_seqs.as_ref().and_then(|v| v.get(i));
            let c = match lookup_tag_coverage(*tag, seq, sample_counts, mismatch) {
                Some(c) if c > 0 => c,
                _ => continue,
            };
            // winner table 目前基于 exact tag hash 构建；mismatch 命中时若 exact hash 不在 winner_map 中，
            // 则保守地视为归属于当前基因组（不扣除），避免过度惩罚真实但带 1 mismatch 的 tag。
            match winner_map.get(tag) {
                Some(winner_entry) if Some(winner_entry.genome_idx) != my_idx => {
                    tags_lost_count += 1; // 该 tag 被分配给了其他基因组
                }
                _ => {
                    covs.push(c); // 归当前基因组所有
                }
            }
        }

        let total_ref_tags = db_entry.tags.len();

        if log_reassign && tags_lost_count > 0 {
            eprintln!(
                "Genome {} in sample {}: owned_tags={}, lost_tags={}, total_ref_tags={}",
                genome_id,
                sample_source,
                covs.len(),
                tags_lost_count,
                total_ref_tags
            );
        }

        // 覆盖度校正后的统计
        let mut result = calculate_statistics(
            covs,
            &db_entry.tag_lengths,
            total_sample_tags,
            total_ref_tags,
            mismatch,
            error_rate,
        );
        result.sample_file = sample_source.to_string();
        result.genome_file = genome_id.to_string();
        result.contig_name = db_entry.sequence_id.clone();

        // 应用profile专用的过滤条件
        if filter_results_for_profile(&result, Some(min_ani), min_shared_tags, min_tags_genome) {
            results.push(result);
        } else if reassign_protection {
            // 过度剥离保护：winner-take-all 在密集多酶库上会把近缘菌株的整个 tag 集判给
            // 赢家，导致输家 post-reassignment ANI/shared_tags 崩塌并被 min-ani 过滤掉
            // （此时 filter_over_reassigned_genomes 根本来不及触发）。触发条件（v2，
            // 两个条件同时满足）：
            //   (a) 近全灭：重分配后存活的命中 tag ≤ REASSIGN_PROTECT_SURVIVAL_FRAC ×
            //       初始共享 tag 数（真阳性湮灭的特征是 2312/2312 全丢；普通菌株簇
            //       输家只丢一部分，不应被保护 —— v1 的 >50% 丢失阈值在此非特异）；
            //   (b) 强初始证据：初始共享 tag 数 ≥ max(REASSIGN_PROTECT_MIN_INITIAL_TAGS,
            //       min_shared_tags)。
            // 满足则以**重分配前**的 ANI/containment/tag 数保留该基因组。
            //
            // 设计取舍：
            // - ANI/tag 数取初始值（检测证据不该被重分配摧毁），但 eff_cov 取
            //   **重分配后**的值 —— 丰度只计入该基因组实际赢得的 tag，避免与赢家
            //   双重计数（重分配只能调丰度，不能杀死检测）。
            // - 防膨胀护栏：保留的基因组必须用其初始统计重新通过 profile 过滤
            //   （min-ani 作用于初始 ANI，min_shared_tags 作用于初始共享数）；
            //   初始 pass 本来已过滤过，这里显式复查以确保不变量。
            if let Some(init) = initial_by_contig.get(db_entry.sequence_id.as_str()) {
                let initial_shared = init.shared_tags;
                let survival_frac = result.shared_tags as f64 / initial_shared.max(1) as f64;
                let min_initial = min_shared_tags.max(REASSIGN_PROTECT_MIN_INITIAL_TAGS);
                if survival_frac <= REASSIGN_PROTECT_SURVIVAL_FRAC
                    && initial_shared >= min_initial
                    && filter_results_for_profile(
                        init,
                        Some(min_ani),
                        min_shared_tags,
                        min_tags_genome,
                    )
                {
                    let mut protected = (*init).clone();
                    protected.sample_file = sample_source.to_string();
                    protected.genome_file = genome_id.to_string();
                    protected.contig_name = db_entry.sequence_id.clone();
                    // 丰度相关字段用重分配后的值（只计实际赢得的 tag）
                    protected.eff_cov = result.eff_cov;
                    protected.eff_lambda = result.eff_lambda;
                    eprintln!(
                        "Reassignment protection: genome {} in sample {} kept only {}/{} initial shared tags ({:.1}% survival) and was stripped below the reporting floor; keeping pre-reassignment ANI {:.2} with post-reassignment eff_cov {:.4} (reassignment-contested)",
                        genome_id, sample_source, result.shared_tags, initial_shared,
                        survival_frac * 100.0, protected.adjusted_ani, protected.eff_cov
                    );
                    results.push(protected);
                }
            }
        }
    }

    results
}

// 使用winner table重新计算统计结果 - 模仿sylph的get_stats函数
//
// `initial_results`：重分配前（初始 pass）的结果，用于过度剥离保护
// （reassign_protection = true 时）。初始结果按 contig_name（= 聚合后的基因组 id）索引。
#[allow(clippy::too_many_arguments)]
fn recalculate_with_winner_table(
    db_entries: &[SyldbEntry],
    sample_entries: &[SylspEntry],
    winner_map: &FxHashMap<Hash, WinnerTableEntry>,
    interner: &GenomeInterner,
    min_ani: f64,
    log_reassign: bool,
    mismatch: usize,
    min_shared_tags: usize,
    min_tags_genome: usize,
    read_error_rate: Option<f64>,
    no_error_correction: bool,
    initial_results: &[QueryResult],
    reassign_protection: bool,
) -> Vec<QueryResult> {
    eprintln!(
        "Recalculating with winner table for {} database entries and {} sample entries",
        db_entries.len(),
        sample_entries.len()
    );

    // 初始结果索引：contig_name（聚合后 = 基因组 id）-> 初始结果。
    // 注意不能用 genome_file 做 key：初始 pass 里 genome_file 被设为 db 路径（所有结果相同）。
    let initial_by_contig: FxHashMap<&str, &QueryResult> = initial_results
        .iter()
        .map(|r| (r.contig_name.as_str(), r))
        .collect();

    // 按样本源分组样本条目
    let mut sample_groups: FxHashMap<String, Vec<&SylspEntry>> = FxHashMap::default();
    for entry in sample_entries {
        sample_groups
            .entry(entry.sample_source.clone())
            .or_default()
            .push(entry);
    }

    let mut all_results = Vec::new();
    
    // 为每个样本源分别计算
    for (sample_source, group_entries) in sample_groups {
        let sample_counts = sample_tag_counts_ref(&group_entries);
        let total_sample_tags = group_entries.len();
        let mut results = recalculate_one_sample_with_winner_table(
            &sample_source,
            &sample_counts,
            total_sample_tags,
            db_entries,
            winner_map,
            interner,
            min_ani,
            log_reassign,
            mismatch,
            min_shared_tags,
            min_tags_genome,
            read_error_rate,
            no_error_correction,
            &initial_by_contig,
            reassign_protection,
        );
        all_results.append(&mut results);
    }

    eprintln!(
        "Recalculation completed: {} results after filtering",
        all_results.len()
    );
    all_results
}

// 过滤过度重新分配的基因组 - 完全模仿sylph的derep_if_reassign_threshold函数
fn filter_over_reassigned_genomes(
    results_old: &[QueryResult],
    results_new: &[QueryResult],
    ani_thresh: f64,
    k: f64
) -> Vec<QueryResult> {
    eprintln!("Filtering over-reassigned genomes: old_results={}, new_results={}, ani_thresh={}, k={}", 
              results_old.len(), results_new.len(), ani_thresh, k);
    
    let ani_thresh = ani_thresh / 100.0;
    
    // 构建genome_id到旧结果的映射
    let mut genome_to_old_result: FxHashMap<String, &QueryResult> = FxHashMap::default();
    for result in results_old.iter() {
        genome_to_old_result.insert(result.genome_file.clone(), result);
    }
    
    let threshold = f64::powf(ani_thresh, k);
    let mut return_vec = Vec::new();
    let mut filtered_count = 0;
    
    for result in results_new.iter() {
        if let Some(old_res) = genome_to_old_result.get(&result.genome_file) {
            let num_tag_reassign = (old_res.shared_tags - result.shared_tags) as f64;
            let reass_thresh = threshold * result.ref_tags as f64;
            
            if num_tag_reassign < reass_thresh {
                return_vec.push(result.clone());
            } else {
                eprintln!("genome {} had num tags reassigned = {}, threshold was {}, removing.", 
                         result.genome_file, num_tag_reassign, reass_thresh);
                filtered_count += 1;
            }
        } else {
            // 如果没有找到旧结果，保留新结果
            return_vec.push(result.clone());
        }
    }
    
    eprintln!("Filtering completed: {} genomes filtered out, {} genomes retained", 
              filtered_count, return_vec.len());
    
    return_vec
}

// 计算 G-score 的函数
// 输入：某个物种的 reads 数 S（所有 2bRAD markers 的总 reads）和测得的 tag 数目 t
// 输出：G-score = sqrt(S * t)
fn calculate_gscore(reads_count: usize, tag_count: usize) -> f64 {
    ((reads_count as f64) * (tag_count as f64)).sqrt()
}

// 基于 G-score 过滤的函数
// 输入：物种列表（包含 S 和 t 信息）以及一个外部指定的阈值 gscore_threshold
// 输出：过滤后的物种列表，只保留 G-score >= gscore_threshold 的物种
fn filter_species_by_gscore(
    species_results: &mut Vec<SpeciesAbundanceResult>,
    gscore_threshold: f64
) -> Vec<SpeciesAbundanceResult> {
    eprintln!("Filtering species by G-score threshold: {:.2}", gscore_threshold);
    
    let initial_count = species_results.len();
    
    // 首先计算每个物种的 G-score
    for species in species_results.iter_mut() {
        species.gscore = calculate_gscore(species.reads_count, species.total_tags);
    }
    
    // 过滤 G-score >= gscore_threshold 的物种
    let filtered_results: Vec<SpeciesAbundanceResult> = species_results
        .iter()
        .filter(|species| species.gscore >= gscore_threshold)
        .cloned()
        .collect();
    
    let filtered_count = filtered_results.len();
    let removed_count = initial_count - filtered_count;
    
    eprintln!("G-score filtering results: {} species retained, {} species removed (threshold: {:.2})", 
              filtered_count, removed_count, gscore_threshold);
    
    if removed_count > 0 {
        eprintln!("Removed species had G-scores below {:.2}", gscore_threshold);
    }
    
    filtered_results
}

// 重新计算丰度 - 模仿sylph的丰度计算逻辑
fn recalculate_abundances_after_reassignment(
    results: &mut [QueryResult],
    sample_entries: &[SylspEntry]
) {
    eprintln!("Recalculating abundances for {} results", results.len());
    
    // 计算总覆盖度
    let total_cov: f64 = results.iter()
        .map(|r| r.eff_cov)
        .sum();
    
    let total_seq_cov: f64 = results.iter()
        .map(|r| r.eff_cov * r.ref_tags as f64)
        .sum();
    
    eprintln!("Total coverage: {:.6}, Total sequence coverage: {:.6}", total_cov, total_seq_cov);
    
    // 重新计算每个结果的丰度
    for result in results.iter_mut() {
        if total_cov > 0.0 {
            result.taxonomic_abundance = result.eff_cov / total_cov * 100.0;
        } else {
            result.taxonomic_abundance = 0.0;
        }
        
        if total_seq_cov > 0.0 {
            result.sequence_abundance = result.eff_cov * result.ref_tags as f64 / total_seq_cov * 100.0;
        } else {
            result.sequence_abundance = 0.0;
        }
    }
    
    eprintln!("Abundance recalculation completed");
}

// 阈值（对齐 sylph 默认值）
const MIN_ANI: f64 = 90.0;                // sylph query 默认 (was 95)
const PROFILE_MIN_ANI: f64 = 95.0;        // sylph profile 默认 (was 97)
pub const MIN_TAGS_FOR_GENOME: usize = 50;    // 基因组最小标签数 (sylph min_number_kmers)
pub const MIN_SHARED_TAGS_FOR_PROFILE: usize = 100; // 覆盖度校正后报告一个基因组所需的最小共享 tag 数；防止极低 tag 数下的虚假优势调用
/// Winner-take-all 过度剥离保护阈值（v2，按 CAMI2 143k 基因组库实测调紧）：
/// 仅当某基因组重分配后**存活**的命中 tag 比例 ≤ 此值（近乎全灭）才考虑保护。
/// v1 用"丢失 >50%"做触发，在密集库上非特异 —— 普通菌株簇的输家也会丢一半 tag，
/// 导致 ~44k 基因组/样本被保护、precision 崩塌（0.69→0.25）。真阳性湮灭的特征是
/// **近全灭**（实测 2312/2312、7326/7327、1555→0、462→1），故改为存活率 ≤5%。
const REASSIGN_PROTECT_SURVIVAL_FRAC: f64 = 0.05;
/// 保护还要求初始共享 tag 数 ≥ 此值（强初始证据）。有效下限为
/// max(REASSIGN_PROTECT_MIN_INITIAL_TAGS, min_shared_tags)。初始证据弱的基因组
/// 即使被全灭也不保护 —— 它们更可能是低质量匹配而非被误杀的真阳性。
const REASSIGN_PROTECT_MIN_INITIAL_TAGS: usize = 100;
// 2bRAD tag 长度，用作 containment->ANI 的指数。BcgI=32；理想情况下应随酶/数据库存储。
const K: f64 = 32.0;

// sylph 覆盖度校正参数（来自 sylph constants.rs / inference.rs）
const SAMPLE_SIZE_CUTOFF: usize = 25;     // 估计 lambda 所需的最小非零覆盖点数
const MEDIAN_ANI_THRESHOLD: f64 = 2.0;    // 中位覆盖 > 此值时朴素 ANI 已足够准确
const MIN_COUNT_CORRECT: f64 = 3.0;       // sylph ratio_lambda 的 min_count_correct 默认值

struct MultiWriter {
    writers: Vec<Box<dyn Write + Send>>,
}

impl MultiWriter {
    fn new() -> Self {
        MultiWriter { writers: Vec::new() }
    }
    fn add_writer(&mut self, writer: Box<dyn Write + Send>) {
        self.writers.push(writer);
    }
}

impl Write for MultiWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        for w in &mut self.writers {
            w.write_all(buf)?;
        }
        Ok(buf.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        for w in &mut self.writers {
            w.flush()?;
        }
        Ok(())
    }
}

pub fn query(args: ContainArgs) -> Result<()> {
    // sketch 模式下 mismatch 需要 tag 序列生成 1-mismatch 邻居，而 k-mer sketch 不存序列。
    if args.sketch_mode && args.mismatch != 0 {
        return Err(anyhow!(
            "--sketch-mode requires --mismatch 0: k-mer sketches (.db/.sp) do not store tag sequences, so 1-mismatch neighbor generation is impossible"
        ));
    }

    // 首先测试文件格式
    let db_files: Vec<_> = args.files.iter().filter(|f| f.ends_with(".db")).collect();

    let sample_files: Vec<_> = args.files.iter().filter(|f| f.ends_with(".sp")).collect();

    if db_files.is_empty() {
        return Err(anyhow!("No .db files found in input files"));
    }

    if sample_files.is_empty() {
        return Err(anyhow!("No .sp files found in input files"));
    }

    // 创建输出写入器
    let writer = Arc::new(Mutex::new(create_multi_writer(&args.out_file_name)?));

    // 打印表头（只打印一次）
    print_header(&writer)?;

    // 并行处理所有数据库文件
    for db_path in db_files {
        eprintln!("Processing database file: {}", db_path);
        
        // 读取数据库文件
        let db_file = File::open(db_path)
            .with_context(|| format!("Failed to open database file: {}", db_path))?;
        let db_reader = BufReader::new(db_file);
        let db_entries: Vec<SyldbEntry> = if args.sketch_mode {
            // sketch 模式：.db 是 Vec<GenomeSketch>，逐条适配成 SyldbEntry
            let genome_sketches: Vec<GenomeSketch> = bincode::deserialize_from(db_reader)
                .with_context(|| format!("Failed to deserialize sketch database file (--sketch-mode expects Vec<GenomeSketch>): {}", db_path))?;
            genome_sketches
                .iter()
                .map(genome_sketch_to_syldb_entry)
                .collect()
        } else {
            bincode::deserialize_from(db_reader)
                .with_context(|| format!("Failed to deserialize database file: {}", db_path))?
        };
        // 关键：按基因组聚合所有 contig（去重 tag），避免逐 contig 碎片化
        let db_entries = aggregate_db_by_genome(db_entries);

        eprintln!("Found {} genomes in database (after aggregating contigs)", db_entries.len());
        let db_enzyme = detect_db_enzyme(&db_entries);
        if let Some(enzyme) = &db_enzyme {
            eprintln!("Database enzyme: {}", enzyme);
        }
        ensure_tag_seqs_for_mismatch(&db_entries, args.mismatch)?;

        // 并行处理所有样本文件
        sample_files
            .par_iter()
            .try_for_each(|sample_path| -> Result<()> {
                eprintln!("Processing sample file: {}", sample_path);

                // sketch 模式：.sp 是 SequencesSketch（或合并的 Vec<SequencesSketch>），
                // 直接借用 kmer_counts，不展开成逐条观测。
                if args.sketch_mode {
                    let sketches = load_sequences_sketches(sample_path)?;
                    for sketch in &sketches {
                        if sketch.kmer_counts.is_empty() {
                            eprintln!(
                                "Warning: Sample sketch {} in {} has no k-mers",
                                sketch.file_name, sample_path
                            );
                            continue;
                        }
                        let sample_enzyme = sketch_enzyme_tag(sketch.c, sketch.k);
                        check_enzyme_compat(
                            db_enzyme.as_deref(),
                            Some(&sample_enzyme),
                            args.allow_enzyme_mismatch,
                            sample_path,
                        )?;

                        let (sample_source, sample_counts, total_sample_tags) =
                            sketch_sample_view(sketch);
                        eprintln!(
                            "Sample {}: {} distinct k-mers, {} total observations",
                            sample_source,
                            sample_counts.len(),
                            total_sample_tags
                        );
                        let error_rate = sample_error_rate(
                            &db_entries,
                            sample_counts,
                            args.mismatch,
                            args.read_error_rate,
                            args.no_error_correction,
                        );

                        // 对每个基因组进行比对
                        for db_entry in &db_entries {
                            let covs: Vec<u32> = db_entry
                                .tags
                                .iter()
                                .enumerate()
                                .filter_map(|(i, t)| {
                                    let seq = db_entry.tag_seqs.as_ref().and_then(|v| v.get(i));
                                    lookup_tag_coverage(*t, seq, sample_counts, args.mismatch)
                                })
                                .collect();

                            // 覆盖度校正后的 ANI 统计
                            let mut result = calculate_statistics(
                                covs,
                                &db_entry.tag_lengths,
                                total_sample_tags,
                                db_entry.tags.len(),
                                args.mismatch,
                                error_rate,
                            );
                            result.sample_file = sample_source.clone();
                            result.genome_file = db_path.to_string();
                            result.contig_name = db_entry.sequence_id.clone();

                            // 应用过滤条件
                            if filter_results(&result, args.minimum_ani, args.min_tags_genome) {
                                print_result(&result, &writer)?;
                            }
                        }
                    }
                    return Ok(());
                }

                let sample_file = File::open(sample_path)
                    .with_context(|| format!("Failed to open sample file: {}", sample_path))?;
                let sample_reader = BufReader::new(sample_file);
                let sample_entries: Vec<SylspEntry> = bincode::deserialize_from(sample_reader)
                    .with_context(|| {
                        format!("Failed to deserialize sample file: {}", sample_path)
                    })?;

                eprintln!("Found {} entries in sample", sample_entries.len());

                // 检查样本数据的有效性
                if sample_entries.is_empty() {
                    eprintln!("Warning: Sample {} has no tags", sample_path);
                    return Ok(());
                }

                // 酶集合一致性检查（DB vs 样本 sketch）
                check_enzyme_compat(
                    db_enzyme.as_deref(),
                    detect_sample_enzyme(&sample_entries).as_deref(),
                    args.allow_enzyme_mismatch,
                    sample_path,
                )?;

                // 统计样本中每个 tag 的覆盖度（重数），用于覆盖度校正
                let sample_counts = sample_tag_counts(&sample_entries);
                let total_sample_tags = sample_entries.len();
                eprintln!("Total tags in sample: {}", total_sample_tags);
                let error_rate = sample_error_rate(
                    &db_entries,
                    &sample_counts,
                    args.mismatch,
                    args.read_error_rate,
                    args.no_error_correction,
                );

                // 对每个基因组进行比对
                for db_entry in &db_entries {
                    // 命中 tag 的覆盖度向量（支持 ≤1 mismatch）
                    let covs: Vec<u32> = db_entry
                        .tags
                        .iter()
                        .enumerate()
                        .filter_map(|(i, t)| {
                            let seq = db_entry.tag_seqs.as_ref().and_then(|v| v.get(i));
                            lookup_tag_coverage(*t, seq, &sample_counts, args.mismatch)
                        })
                        .collect();

                    // 覆盖度校正后的 ANI 统计
                    let mut result = calculate_statistics(
                        covs,
                        &db_entry.tag_lengths,
                        total_sample_tags,
                        db_entry.tags.len(),
                        args.mismatch,
                        error_rate,
                    );
                    result.sample_file = sample_path.to_string();
                    result.genome_file = db_path.to_string();
                    result.contig_name = db_entry.sequence_id.clone();

                    // 应用过滤条件
                    if filter_results(&result, args.minimum_ani, args.min_tags_genome) {
                        print_result(&result, &writer)?;
                    }
                }
                Ok(())
            })?;
    }

    Ok(())
}

fn create_multi_writer(out_file_name: &Option<String>) -> Result<Box<dyn Write + Send>> {
    let mut mw = MultiWriter::new();
    mw.add_writer(Box::new(BufWriter::new(std::io::stdout())));
    if let Some(path) = out_file_name {
        let file = File::create(path)
            .with_context(|| format!("Failed to create output file: {}", path))?;
        mw.add_writer(Box::new(BufWriter::new(file)));
    }
    Ok(Box::new(mw))
}

fn print_header(writer: &Arc<Mutex<Box<dyn Write + Send>>>) -> Result<()> {
    let mut writer = writer.lock().unwrap();
    writeln!(writer, "{:<20} {:<20} {:<10} {:<10} {:<15} {:<15} {:<10} {:<10} {:<10} {:<15} {:<10} {:<10}",
        "Sample_file", "Genome_file", "ANI(%)", "Eff_cov", "ANI_5-95%", "Eff_lambda", "Lambda_5-95%", "Median_cov", "Mean_cov", "Containment", "Naive_ANI", "Contig_name")?;
    writeln!(writer, "{:-<150}", "")?;
    Ok(())
}

fn print_result(result: &QueryResult, writer: &Arc<Mutex<Box<dyn Write + Send>>>) -> Result<()> {
    let mut writer = writer.lock().unwrap();
    writeln!(writer, "{:<20} {:<20} {:<10.2} {:<10.3} {:<7.2}-{:<7.2} {:<10.3} {:<7.2}-{:<7.2} {:<10.3} {:<10.3} {:<7} {:<10.2} {:<10}",
        result.sample_file,
        result.genome_file,
        result.adjusted_ani,
        result.eff_cov,
        result.ani_percentile.0,
        result.ani_percentile.1,
        result.eff_lambda,
        result.lambda_percentile.0,
        result.lambda_percentile.1,
        result.median_cov,
        result.mean_cov_geq1,
        result.containment_ind,
        result.naive_ani,
        result.contig_name
    )?;
    Ok(())
}

fn mean_u32(v: &[u32]) -> f64 {
    if v.is_empty() { 0.0 } else { v.iter().map(|&x| x as f64).sum::<f64>() / v.len() as f64 }
}

/// sylph 的 ratio 估计器：从覆盖度直方图估计泊松均值 lambda。
fn ratio_lambda(full_covs: &[u32], min_count_correct: f64) -> Option<f64> {
    let mut num_zero = 0usize;
    let mut count_map: FxHashMap<usize, usize> = FxHashMap::default();
    for &x in full_covs {
        if x == 0 { num_zero += 1; }
        else { *count_map.entry(x as usize).or_insert(0) += 1; }
    }
    if count_map.len() <= 1 { return None; }
    if full_covs.len() - num_zero < SAMPLE_SIZE_CUTOFF { return None; }
    let mut sort_vec: Vec<(usize, usize)> = count_map.iter().map(|(&k, &v)| (v, k)).collect();
    sort_vec.sort_by(|a, b| b.0.cmp(&a.0));
    let most_ind = sort_vec[0].1;
    let count_p1 = *count_map.get(&(most_ind + 1))? as f64;
    let count = count_map[&most_ind] as f64;
    if count_p1 < min_count_correct || count < min_count_correct { return None; }
    Some(count_p1 / count * (most_ind + 1) as f64)
}

/// sylph 的覆盖度校正 ANI：用零截断泊松的检出概率 (1 - e^-lambda) 还原因低覆盖
/// 而漏掉的 tag，再映射回 ANI。这是 meta2bseek 之前缺失、导致大量假阴性的核心。
fn ani_from_lambda(lambda: f64, k: f64, full_covs: &[u32]) -> Option<f64> {
    if full_covs.is_empty() { return None; }
    let contain_count = full_covs.iter().filter(|&&x| x != 0).count() as f64;
    let denom = 1.0 - (-lambda).exp();
    if denom <= 0.0 { return None; }
    let adj_index = contain_count / denom / full_covs.len() as f64;
    let ani = adj_index.powf(1.0 / k);
    if ani.is_nan() || ani < 0.0 { None } else { Some(ani) }
}

fn empty_result(shared: usize, total_ref_tags: usize, query_tags: usize) -> QueryResult {
    QueryResult {
        sample_file: String::new(),
        genome_file: String::new(),
        contig_name: String::new(),
        adjusted_ani: 0.0,
        eff_cov: 0.0,
        ani_percentile: (0.0, 0.0),
        eff_lambda: 0.0,
        lambda_percentile: (0.0, 0.0),
        median_cov: 0.0,
        mean_cov_geq1: 0.0,
        containment_ind: format!("{}/{}", shared, total_ref_tags),
        naive_ani: 0.0,
        ref_tags: total_ref_tags,
        shared_tags: shared,
        query_tags,
        taxonomic_abundance: 0.0,
        sequence_abundance: 0.0,
    }
}

/// 多酶联合 ANI：按 tag 长度分区，信息加权最小二乘估计。
/// `classes` 映射 length -> (N_l, D_l)；`lambda` 为覆盖度校正参数，None 表示不做校正。
/// `e` 为 per-base read error rate（仅 mismatch=1 的 error-aware 反推使用；e=0 即旧行为）。
fn joint_ani_from_classes(classes: &FxHashMap<u8, (usize, usize)>, lambda: Option<f64>, mismatch: usize, e: f64) -> Option<f64> {
    if mismatch == 0 {
        // exact matching：c_adj = a^ℓ，对 log 空间做加权最小二乘。
        let mut num = 0.0;
        let mut den = 0.0;
        for (&len, &(n_l, d_l)) in classes {
            if n_l == 0 || d_l == 0 {
                continue;
            }
            let c = d_l as f64 / n_l as f64;
            let c_adj = if let Some(lam) = lambda {
                let p = 1.0 - (-lam).exp();
                if p <= 0.0 {
                    continue;
                }
                (c / p).min(1.0)
            } else {
                c
            };
            let l = len as f64;
            num += n_l as f64 * l * c_adj.ln();
            den += n_l as f64 * l * l;
        }
        if den == 0.0 {
            return None;
        }
        Some((num / den).exp())
    } else {
        // ≤1 mismatch：先除以覆盖度校正因子得到 P_detect(q,ℓ)，再逐长度类反解 q 并除以
        // (1-e) 恢复真实 ANI a（error-aware；e=0 时与旧公式完全一致），最后加权平均。
        let mut weighted_a = 0.0;
        let mut total_weight = 0.0;
        for (&len, &(n_l, d_l)) in classes {
            if n_l == 0 || d_l == 0 {
                continue;
            }
            let c = d_l as f64 / n_l as f64;
            let c_adj = if let Some(lam) = lambda {
                let p = 1.0 - (-lam).exp();
                if p <= 0.0 {
                    continue;
                }
                (c / p).min(1.0)
            } else {
                c
            };
            let a_est = ani_from_containment_one_mismatch_err(c_adj, len as usize, e);
            let weight = n_l as f64 * len as f64;
            weighted_a += weight * a_est;
            total_weight += weight;
        }
        if total_weight == 0.0 {
            return None;
        }
        Some((weighted_a / total_weight).min(1.0))
    }
}

/// 计算单个基因组的朴素 + 覆盖度校正 ANI（sylph get_stats 风格）。
/// `covs` = 该基因组在样本中命中的 tag 的覆盖度（重数），长度 = 共享 tag 数；
/// `tag_lengths` = 该基因组每个参考 tag 的长度（与完整 tag 集一一对应）；
/// `query_tags` = 样本中总 tag 数；
/// `total_ref_tags` = 该基因组去重后的 tag 总数。
fn calculate_statistics(
    mut covs: Vec<u32>,
    tag_lengths: &[u8],
    query_tags: usize,
    total_ref_tags: usize,
    mismatch: usize,
    error_rate: f64,
) -> QueryResult {
    let shared = covs.len();
    if shared == 0 || total_ref_tags == 0 || tag_lengths.len() != total_ref_tags {
        return empty_result(shared, total_ref_tags, query_tags);
    }

    let containment = shared as f64 / total_ref_tags as f64;

    covs.sort_unstable();
    let median_cov = covs[covs.len() / 2] as f64;
    let mean_geq1 = mean_u32(&covs);

    // 完整覆盖向量：命中 tag 的覆盖度 + 漏掉的 tag 记为 0
    let mut full_covs = covs;
    full_covs.extend(std::iter::repeat(0u32).take(total_ref_tags - shared));

    // 按长度统计 (N_l, D_l)
    let mut length_classes: FxHashMap<u8, (usize, usize)> = FxHashMap::default();
    for (i, &len) in tag_lengths.iter().enumerate() {
        let entry = length_classes.entry(len).or_insert((0, 0));
        entry.0 += 1;
        if full_covs[i] > 0 {
            entry.1 += 1;
        }
    }

    // 朴素 ANI：不做覆盖度校正
    let naive_ani = joint_ani_from_classes(&length_classes, None, mismatch, error_rate)
        .unwrap_or_else(|| {
            if mismatch == 1 {
                // 单长度类时退到主导长度；无长度信息时用 K=32 近似。
                let dominant_len = length_classes.iter()
                    .max_by_key(|(_, (n, _))| *n)
                    .map(|(&l, _)| l as usize)
                    .unwrap_or(K as usize);
                ani_from_containment_one_mismatch_err(containment, dominant_len, error_rate)
            } else {
                containment.powf(1.0 / K)
            }
        });

    // 是否需要覆盖度校正
    let (adjusted_ani, eff_cov) = if median_cov > MEDIAN_ANI_THRESHOLD {
        // 覆盖足够高，朴素 containment 已能反映一致性
        (naive_ani, median_cov)
    } else {
        match ratio_lambda(&full_covs, MIN_COUNT_CORRECT) {
            Some(lam) => (
                joint_ani_from_classes(&length_classes, Some(lam), mismatch, error_rate).unwrap_or(naive_ani),
                lam,
            ),
            None => (naive_ani, mean_geq1), // 点数不足，无法估计 lambda：回退到朴素 ANI
        }
    };
    let adjusted_ani = adjusted_ani.min(1.0); // 截断到 100%

    let mut result = empty_result(shared, total_ref_tags, query_tags);
    result.naive_ani = naive_ani * 100.0;
    result.adjusted_ani = adjusted_ani * 100.0;
    result.eff_cov = eff_cov;
    result.eff_lambda = eff_cov;
    result.median_cov = median_cov;
    result.mean_cov_geq1 = mean_geq1;

    // 简单置信带（仅用于显示）
    let unc = 1.0 + (1.0 - containment) * 1.5;
    result.ani_percentile = ((result.adjusted_ani - unc).max(0.0), (result.adjusted_ani + unc).min(100.0));
    result.lambda_percentile = ((eff_cov - 0.1).max(0.0), eff_cov + 0.1);

    result
}

/// 把每条 contig 的 SyldbEntry 按基因组（genome_source 推导出的 id）合并成一个
/// SyldbEntry，tag 去重。这样 query/profile 都按"整基因组"而非"逐 contig"统计，
/// 这是降低假阴性最关键的一步。
/// 从 DB 条目中推断构建时使用的酶（取非空众数）。
fn detect_db_enzyme(entries: &[SyldbEntry]) -> Option<String> {
    let mut counts: FxHashMap<String, usize> = FxHashMap::default();
    for e in entries {
        if !e.enzyme.is_empty() {
            *counts.entry(e.enzyme.clone()).or_insert(0) += 1;
        }
    }
    counts.into_iter().max_by_key(|(_, c)| *c).map(|(e, _)| e)
}

/// 从样本条目（.sylsp）中推断提取时使用的酶（取非空众数）。旧格式无此字段 → None。
fn detect_sample_enzyme(entries: &[SylspEntry]) -> Option<String> {
    let mut counts: FxHashMap<String, usize> = FxHashMap::default();
    for e in entries {
        if !e.enzyme.is_empty() {
            *counts.entry(e.enzyme.clone()).or_insert(0) += 1;
        }
    }
    counts.into_iter().max_by_key(|(_, c)| *c).map(|(e, _)| e)
}

/// 把逗号分隔的酶字符串规范化为有序集合（去空白、去重、排序）。
fn enzyme_set(s: &str) -> std::collections::BTreeSet<String> {
    s.split(',')
        .map(|x| x.trim())
        .filter(|x| !x.is_empty())
        .map(String::from)
        .collect()
}

/// DB 与样本的酶集合一致性检查。
/// - 集合相同（或任一侧缺失酶信息，如旧格式文件）→ 直接通过。
/// - 集合不同 → 打印醒目的 stderr 警告，列出两侧酶集合，并提示 DB 中缺失于样本的酶
///   会导致 containment 被压低（这些 tag 在样本中根本不存在）。
/// - 若 DB 含有样本完全缺失的酶：默认硬错误；`--allow-enzyme-mismatch` 可降级为警告继续。
///   （样本多出 DB 没有的酶只是浪费样本 tag，不硬错。）
fn check_enzyme_compat(
    db_enzyme: Option<&str>,
    sample_enzyme: Option<&str>,
    allow_mismatch: bool,
    context: &str,
) -> Result<()> {
    let (db_str, sp_str) = match (db_enzyme, sample_enzyme) {
        (Some(d), Some(s)) if !d.is_empty() && !s.is_empty() => (d, s),
        _ => return Ok(()), // 任一侧无酶信息（旧格式），无法检查
    };
    let db_set = enzyme_set(db_str);
    let sp_set = enzyme_set(sp_str);
    if db_set == sp_set {
        return Ok(());
    }
    let missing: Vec<String> = db_set.difference(&sp_set).cloned().collect();
    let extra: Vec<String> = sp_set.difference(&db_set).cloned().collect();
    eprintln!("================ ENZYME MISMATCH ================");
    eprintln!("{}: database enzymes = {{{}}}, sample enzymes = {{{}}}",
              context,
              db_set.iter().cloned().collect::<Vec<_>>().join(","),
              sp_set.iter().cloned().collect::<Vec<_>>().join(","));
    if !missing.is_empty() {
        eprintln!("DB enzymes MISSING from sample sketch: {{{}}}. \
                  Containment/ANI for tags of these enzymes will be deflated (they can never match). \
                  Re-extract the sample with the same enzyme set as the database.",
                  missing.join(","));
    }
    if !extra.is_empty() {
        eprintln!("Sample enzymes not in database: {{{}}} (those sample tags are simply unused).",
                  extra.join(","));
    }
    eprintln!("=================================================");
    if !missing.is_empty() && !allow_mismatch {
        return Err(anyhow!(
            "enzyme mismatch in {}: database requires {{{}}} but the sample sketch was extracted without them. \
             Re-extract the sample with matching enzymes, or pass --allow-enzyme-mismatch to proceed anyway",
            context, missing.join(",")
        ));
    }
    if !missing.is_empty() {
        eprintln!("--allow-enzyme-mismatch set; continuing despite missing sample enzymes {{{}}}", missing.join(","));
    }
    Ok(())
}

/// mismatch 模式需要数据库中存储的 tag 序列来生成 1-mismatch 邻居；
/// 用 --no-tag-seqs 构建的（或更旧版本的）数据库不含序列，应明确报错而非静默退化为 exact match。
fn ensure_tag_seqs_for_mismatch(entries: &[SyldbEntry], mismatch: usize) -> Result<()> {
    if mismatch > 0 && entries.iter().any(|e| e.tag_seqs.is_none()) {
        return Err(anyhow!(
            "mismatch mode requires a database built with tag sequences; this database was built with --no-tag-seqs (or predates tag sequence storage). Rebuild the database without --no-tag-seqs or use --mismatch 0"
        ));
    }
    Ok(())
}

// ---- sketch-mode adapters: run k-mer sketches through the 2bRAD profiling engine ----

/// Self-describing "enzyme" string for sketch-mode inputs. Both the database and the
/// sample side derive it with the same rule, so the existing `check_enzyme_compat`
/// passes natively when c and k match (and hard-errors when they do not).
fn sketch_enzyme_tag(c: usize, k: usize) -> String {
    format!("kmer_c{}_k{}", c, k)
}

/// Adapt a k-mer `GenomeSketch` into a `SyldbEntry` so the entire profiling engine
/// (per-genome aggregation, containment/ANI statistics, winner-take-all reassignment)
/// is reused unchanged. `aggregate_db_by_genome` downstream deduplicates `tags`.
/// Sketches store no tag sequences, so `tag_seqs` is None and only --mismatch 0 is valid.
fn genome_sketch_to_syldb_entry(sketch: &GenomeSketch) -> SyldbEntry {
    SyldbEntry {
        sequence_id: sketch.file_name.clone(),
        tags: sketch.genome_kmers.clone(),
        tag_lengths: vec![sketch.k as u8; sketch.genome_kmers.len()],
        genome_source: sketch.file_name.clone(),
        tag_uniqueness: None,
        tag_seqs: None,
        enzyme: sketch_enzyme_tag(sketch.c, sketch.k),
    }
}

/// Load a sketch-mode sample file: either a single `SequencesSketch` or a merged
/// `Vec<SequencesSketch>` (same two-form deserialization as `view`).
fn load_sequences_sketches(sample_path: &str) -> Result<Vec<SequencesSketch>> {
    let file = File::open(sample_path)
        .with_context(|| format!("Failed to open sample file: {}", sample_path))?;
    let reader = BufReader::with_capacity(100_000_000, file);
    if let Ok(sketch) = bincode::deserialize_from::<_, SequencesSketch>(reader) {
        return Ok(vec![sketch]);
    }
    // First attempt consumed the reader; reopen before trying the merged form.
    let file = File::open(sample_path)
        .with_context(|| format!("Failed to open sample file: {}", sample_path))?;
    let reader = BufReader::with_capacity(100_000_000, file);
    let sketches: Vec<SequencesSketch> = bincode::deserialize_from(reader)
        .with_context(|| format!(
            "Failed to deserialize sample sketch file (--sketch-mode expects SequencesSketch or Vec<SequencesSketch>): {}",
            sample_path
        ))?;
    Ok(sketches)
}

/// Per-sample view of a `SequencesSketch`: (sample source, k-mer counts, total observations).
/// The counts map is borrowed directly — it is NEVER expanded into per-observation
/// entries (real samples hold tens of millions of observations; expanding them would
/// explode memory).
fn sketch_sample_view(sketch: &SequencesSketch) -> (String, &FxHashMap<Hash, u32>, usize) {
    let source = sketch
        .sample_name
        .clone()
        .unwrap_or_else(|| sketch.file_name.clone());
    let total = sketch.kmer_counts.values().map(|&c| c as usize).sum();
    (source, &sketch.kmer_counts, total)
}

fn aggregate_db_by_genome(entries: Vec<SyldbEntry>) -> Vec<SyldbEntry> {
    let mut map: FxHashMap<String, SyldbEntry> = FxHashMap::default();
    for e in entries {
        let gid = extract_genome_id_from_path(&e.genome_source).to_string();
        let agg = map.entry(gid.clone()).or_insert_with(|| SyldbEntry {
            sequence_id: gid.clone(),
            tags: Vec::new(),
            tag_lengths: Vec::new(),
            genome_source: e.genome_source.clone(),
            tag_uniqueness: None,
            tag_seqs: None,
            enzyme: e.enzyme.clone(),
        });
        agg.tags.extend(e.tags);
        agg.tag_lengths.extend(e.tag_lengths);
        if let Some(e_seqs) = e.tag_seqs.as_ref() {
            agg.tag_seqs.get_or_insert_with(Vec::new).extend(e_seqs.iter().cloned());
        }
    }
    let mut out: Vec<SyldbEntry> = map.into_values().collect();
    for e in &mut out {
        // 按 (hash, length, seq) 去重，保证多酶模式下同一 hash 的不同长度 tag 不被错误合并
        let n = e.tags.len();
        let seqs = e.tag_seqs.as_ref().map(|s| s.as_slice());
        let mut triples: Vec<(Hash, u8, TagHash)> = Vec::with_capacity(n);
        for i in 0..n {
            let seq = seqs.map(|s| s.get(i).cloned()).flatten().unwrap_or_default();
            triples.push((e.tags[i], e.tag_lengths[i], seq));
        }
        triples.sort_unstable();
        triples.dedup();
        e.tags = triples.iter().map(|(h, _, _)| *h).collect();
        e.tag_lengths = triples.iter().map(|(_, l, _)| *l).collect();
        if e.tag_seqs.is_some() {
            e.tag_seqs = Some(triples.iter().map(|(_, _, s)| s.clone()).collect());
        }
    }
    out
}

/// 统计样本中每个 tag 的出现次数（覆盖度/重数）。之前的实现把重数压成 1，
/// 导致无法估计覆盖度 lambda。
fn sample_tag_counts(sample_entries: &[SylspEntry]) -> FxHashMap<Hash, u32> {
    let mut counts: FxHashMap<Hash, u32> = FxHashMap::default();
    for e in sample_entries {
        *counts.entry(e.tag).or_insert(0) += 1;
    }
    counts
}

fn sample_tag_counts_ref(sample_entries: &[&SylspEntry]) -> FxHashMap<Hash, u32> {
    let mut counts: FxHashMap<Hash, u32> = FxHashMap::default();
    for e in sample_entries {
        *counts.entry(e.tag).or_insert(0) += 1;
    }
    counts
}

/// 查询一个参考 tag 在样本中的覆盖度。`mismatch=0` 时精确匹配；`mismatch=1` 时额外检查
/// 该参考 tag 的所有 canonical 1-mismatch 邻居。优先返回 exact match 的计数，若不存在则
/// 返回第一个命中的 neighbor 计数（避免同一参考 tag 被多个样本 tag 重复计数）。
fn lookup_tag_coverage(
    tag_hash: Hash,
    tag_seq: Option<&TagHash>,
    sample_counts: &FxHashMap<Hash, u32>,
    mismatch: usize,
) -> Option<u32> {
    if let Some(&c) = sample_counts.get(&tag_hash) {
        if c > 0 {
            return Some(c);
        }
    }
    if mismatch == 0 {
        return None;
    }
    let seq = tag_seq?;
    for neighbor_hash in one_mismatch_canonical_hashes(seq) {
        if let Some(&c) = sample_counts.get(&neighbor_hash) {
            if c > 0 {
                return Some(c);
            }
        }
    }
    None
}

// ---- read error rate 估计（--mismatch 1 的 error-aware ANI 反推） ----

/// e 估计要求候选基因组至少有这么多个 exact-match tag，否则噪声太大。
const MIN_TAGS_FOR_ERROR_EST: usize = 100;
/// 用于 e 估计的"近同一基因组"阈值：mm0 模型 ANI ≥ 此值的基因组才用于估计。
/// 保守取 0.98：更趋异的基因组会把真实分歧 a 误当成错误 e（a/e 混淆）。
const ERROR_EST_MIN_ANI: f64 = 0.98;
/// 没有基因组达到 ERROR_EST_MIN_ANI 时，退而使用的 top-N 候选（按 exact containment）。
const ERROR_EST_TOP_N: usize = 5;
/// 估计出的 e 上限（超出视为估计失败迹象，截断并告警）。
const MAX_READ_ERROR_RATE: f64 = 0.5;

/// 由按长度类汇总的 (exact 命中 read 数 E_l, mm1 命中 read 数 M_l) 估计 per-base read error rate。
///
/// 模型：对近同一基因组（a≈1），每条 read 上的同源 tag 以 per-base 概率 q = a(1-e) ≈ 1-e
/// 被观察到，故按 read 计数 E[mm1]/E[exact] = ℓ(1-q)/q，即
///   q_l = ℓ·E_l / (ℓ·E_l + M_l)，e_l = M_l / (ℓ·E_l + M_l)。
/// 按 read（覆盖度重数）计数而非 distinct tag：每条 read 是独立抽样，比例不受覆盖度
/// 影响；若按 tag 存在性计数，多覆盖 tag 只要有一条无错 read 就记为 exact，会系统性
/// 低估 e。跨长度类以 exact read 数加权合并为全局 e。
/// 返回 None 表示无法估计（调用方回退 e=0）。
fn pooled_error_rate(classes: &FxHashMap<u8, (u64, u64)>) -> Option<f64> {
    let mut weighted_e = 0.0;
    let mut total_exact = 0u64;
    let mut total_mm1 = 0u64;
    for (&len, &(exact, mm1)) in classes {
        total_exact += exact;
        total_mm1 += mm1;
        if exact == 0 {
            continue;
        }
        let l = len as f64;
        let e_l = mm1 as f64 / (l * exact as f64 + mm1 as f64);
        weighted_e += exact as f64 * e_l;
    }
    if total_exact < MIN_TAGS_FOR_ERROR_EST as u64 {
        return None; // tag 太少，估计不可靠
    }
    if total_mm1 == 0 {
        return Some(0.0); // 没有 mm1-only 命中：e≈0（也可能全基因组 exact，属正常）
    }
    Some(weighted_e / total_exact as f64)
}

/// 从 top 候选基因组估计样本的 per-base read error rate e（--mismatch 1 专用）。
///
/// 候选选择：先取 mm0 模型 ANI ≥ ERROR_EST_MIN_ANI 的基因组（保守，减少 a/e 混淆）；
/// 若无，退到 exact containment 最高的前 ERROR_EST_TOP_N 个。对这些基因组统计每个参考
/// tag 是 exact 命中还是 mm1-only 命中，交给 `pooled_error_rate`。
/// 估计不可行（无候选 / tag 太少）时返回 0.0 并打日志 —— 即退化为旧的未校正行为。
fn estimate_read_error_rate(db_entries: &[SyldbEntry], sample_counts: &FxHashMap<Hash, u32>) -> f64 {
    // 第一遍：计算每个基因组的 exact containment 与 mm0-ANI（按主导长度类）。
    let mut scored: Vec<(f64, &SyldbEntry)> = Vec::new(); // (exact containment, entry)
    for db_entry in db_entries {
        if db_entry.tags.is_empty() || db_entry.tag_seqs.is_none() {
            continue;
        }
        let n_exact = db_entry.tags.iter()
            .filter(|t| sample_counts.get(*t).map_or(false, |&c| c > 0))
            .count();
        let c0 = n_exact as f64 / db_entry.tags.len() as f64;
        if c0 > 0.0 {
            scored.push((c0, db_entry));
        }
    }
    if scored.is_empty() {
        eprintln!("Read error rate: no genomes with exact matches; falling back to e=0 (uncorrected mm1 ANI)");
        return 0.0;
    }

    // 主导长度类（按 tag 数）用于 mm0-ANI 粗估。
    let dominant_len = |entry: &SyldbEntry| -> usize {
        let mut by_len: FxHashMap<u8, usize> = FxHashMap::default();
        for &l in &entry.tag_lengths {
            *by_len.entry(l).or_insert(0) += 1;
        }
        by_len.into_iter().max_by_key(|(_, n)| *n).map(|(l, _)| l as usize).unwrap_or(K as usize)
    };

    let mut selected: Vec<&SyldbEntry> = scored.iter()
        .filter(|(c0, entry)| c0.powf(1.0 / dominant_len(entry) as f64) >= ERROR_EST_MIN_ANI)
        .map(|(_, entry)| *entry)
        .collect();
    if selected.is_empty() {
        // 没有近同一基因组：退到 exact containment 最高的 top-N，并在日志中提示
        // 此时 e 可能与真实分歧 a 混淆（趋异群落中应谨慎解释）。
        scored.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
        selected = scored.iter().take(ERROR_EST_TOP_N).map(|(_, entry)| *entry).collect();
        eprintln!("Read error rate: no genome with mm0-ANI >= {:.0}%; using top {} genomes by exact containment (e may be confounded with true divergence)",
                  ERROR_EST_MIN_ANI * 100.0, selected.len());
    }

    // 第二遍：对候选基因组逐 tag 统计 exact / mm1 命中的 read 数（按长度类）。
    // 一个样本 tag 若是多个候选 ref tag 的邻居会被重复计入 M，这种情况罕见，对比例影响可忽略。
    let mut classes: FxHashMap<u8, (u64, u64)> = FxHashMap::default();
    for entry in &selected {
        let tag_seqs = entry.tag_seqs.as_ref().unwrap();
        for (i, tag) in entry.tags.iter().enumerate() {
            let len = entry.tag_lengths.get(i).copied().unwrap_or(K as u8);
            let c_exact = sample_counts.get(tag).copied().unwrap_or(0) as u64;
            let mut c_mm1 = 0u64;
            if let Some(seq) = tag_seqs.get(i) {
                for h in one_mismatch_canonical_hashes(seq) {
                    if let Some(&c) = sample_counts.get(&h) {
                        c_mm1 += c as u64;
                    }
                }
            }
            let e = classes.entry(len).or_insert((0, 0));
            e.0 += c_exact;
            e.1 += c_mm1;
        }
    }

    match pooled_error_rate(&classes) {
        Some(e) if e > MAX_READ_ERROR_RATE => {
            eprintln!("Read error rate: estimated e={:.4} exceeds {:.2}; clamping (estimate likely unreliable)",
                      e, MAX_READ_ERROR_RATE);
            MAX_READ_ERROR_RATE
        }
        Some(e) => {
            eprintln!("Read error rate: estimated e={:.4} from {} candidate genome(s)", e, selected.len());
            e
        }
        None => {
            eprintln!("Read error rate: too few matched tags to estimate; falling back to e=0 (uncorrected mm1 ANI)");
            0.0
        }
    }
}

/// 决定一次比对实际使用的 read error rate：
/// `--no-error-correction` → 0（旧行为）；`--read-error-rate` → 固定值；
/// 否则（仅 mismatch=1）用 `estimate_read_error_rate` 的数据驱动估计。mismatch=0 恒为 0。
fn sample_error_rate(
    db_entries: &[SyldbEntry],
    sample_counts: &FxHashMap<Hash, u32>,
    mismatch: usize,
    override_e: Option<f64>,
    no_error_correction: bool,
) -> f64 {
    if mismatch == 0 || no_error_correction {
        return 0.0;
    }
    if let Some(e) = override_e {
        eprintln!("Read error rate: using fixed e={:.4} from --read-error-rate", e);
        return e.clamp(0.0, MAX_READ_ERROR_RATE);
    }
    estimate_read_error_rate(db_entries, sample_counts)
}

fn filter_results(result: &QueryResult, min_ani: Option<f64>, min_tags_genome: usize) -> bool {
    // 只按 sylph 的方式过滤：有命中、基因组足够大、且（覆盖度校正后的）ANI 达标。
    if result.shared_tags == 0 {
        return false;
    }
    if result.ref_tags < min_tags_genome {
        return false;
    }
    let effective_min_ani = min_ani.unwrap_or(MIN_ANI);
    result.adjusted_ani >= effective_min_ani
}

// profile 专用过滤（阈值更严，但同样依赖覆盖度校正后的 ANI）
fn filter_results_for_profile(result: &QueryResult, min_ani: Option<f64>, min_shared_tags: usize, min_tags_genome: usize) -> bool {
    if result.shared_tags == 0 {
        return false;
    }
    if result.ref_tags < min_tags_genome {
        return false;
    }
    if result.shared_tags < min_shared_tags {
        return false;
    }
    let effective_min_ani = min_ani.unwrap_or(PROFILE_MIN_ANI);
    result.adjusted_ani >= effective_min_ani
}

// Per-sample comparison against the cached, genome-aggregated database — the core of
// query_single_file_with_cached_db. Shared by the 2bRAD path (counts built from
// SylspEntry lists via sample_tag_counts_ref) and the sketch path (counts borrowed
// directly from SequencesSketch::kmer_counts).
#[allow(clippy::too_many_arguments)]
fn query_one_sample_with_cached_db(
    sample_source: &str,
    sample_counts: &FxHashMap<Hash, u32>,
    total_sample_tags: usize,
    db_path: &str,
    cached_db_entries: &[SyldbEntry],
    min_ani: f64,
    mismatch: usize,
    min_shared_tags: usize,
    min_tags_genome: usize,
    read_error_rate: Option<f64>,
    no_error_correction: bool,
) -> Vec<QueryResult> {
    let error_rate = sample_error_rate(
        cached_db_entries,
        sample_counts,
        mismatch,
        read_error_rate,
        no_error_correction,
    );

    // 并行处理每个基因组（已按基因组聚合）进行比对
    cached_db_entries
        .par_iter()
        .filter_map(|db_entry| {
            // 基因组太小/太碎则跳过（sylph min_number_kmers）
            if db_entry.tags.len() < min_tags_genome {
                return None;
            }

            let covs: Vec<u32> = db_entry
                .tags
                .iter()
                .enumerate()
                .filter_map(|(i, t)| {
                    let seq = db_entry.tag_seqs.as_ref().and_then(|v| v.get(i));
                    lookup_tag_coverage(*t, seq, sample_counts, mismatch)
                })
                .collect();

            let mut result = calculate_statistics(
                covs,
                &db_entry.tag_lengths,
                total_sample_tags,
                db_entry.tags.len(),
                mismatch,
                error_rate,
            );
            result.sample_file = sample_source.to_string();
            result.genome_file = db_path.to_string();
            result.contig_name = db_entry.sequence_id.clone();

            // 应用profile专用的过滤条件
            if filter_results_for_profile(&result, Some(min_ani), min_shared_tags, min_tags_genome)
            {
                Some(result)
            } else {
                None
            }
        })
        .collect()
}

// 内部函数：使用缓存的数据库数据进行查询 - 优化大文件读取
fn query_single_file_with_cached_db(
    sample_path: &str,
    db_path: &str,
    cached_db_entries: &[SyldbEntry],
    cached_sample_entries: &FxHashMap<String, Vec<SylspEntry>>,
    min_ani: f64,
    mismatch: usize,
    min_shared_tags: usize,
    min_tags_genome: usize,
    read_error_rate: Option<f64>,
    no_error_correction: bool,
) -> Result<Vec<QueryResult>> {
    eprintln!("Processing sample file with cached database: {}", sample_path);
    
    // 从缓存中获取样本数据
    let sample_entries = cached_sample_entries.get(sample_path)
        .ok_or_else(|| anyhow!("Sample file not found in cache: {}", sample_path))?;

    eprintln!("Found {} entries in sample", sample_entries.len());

    // 检查样本数据的有效性
    if sample_entries.is_empty() {
        eprintln!("Warning: Sample {} has no tags", sample_path);
        return Ok(Vec::new());
    }

    // 按样本源分组 - 这是关键：处理合并文件中的多个样本
    let mut sample_groups: FxHashMap<String, Vec<&SylspEntry>> = FxHashMap::default();
    for entry in sample_entries {
        sample_groups
            .entry(entry.sample_source.clone())
            .or_default()
            .push(entry);
    }

    eprintln!(
        "Found {} different sample sources in file: {:?}",
        sample_groups.len(),
        sample_groups.keys().collect::<Vec<_>>()
    );

    // 并行处理每个样本组，然后合并结果
    let mut all_results: Vec<QueryResult> = sample_groups
        .par_iter()
        .flat_map(|(sample_source, entries)| {
            eprintln!(
                "Processing sample source: {} with {} entries",
                sample_source,
                entries.len()
            );

            // 统计样本中每个 tag 的覆盖度（重数）
            let sample_counts = sample_tag_counts_ref(entries);
            query_one_sample_with_cached_db(
                sample_source,
                &sample_counts,
                entries.len(),
                db_path,
                cached_db_entries,
                min_ani,
                mismatch,
                min_shared_tags,
                min_tags_genome,
                read_error_rate,
                no_error_correction,
            )
        })
        .collect();

    // 按ANI排序（参考sylph的排序机制）
    all_results.sort_by(|a, b| b.adjusted_ani.partial_cmp(&a.adjusted_ani).unwrap());

    eprintln!("Generated {} results for file {}", all_results.len(), sample_path);
    Ok(all_results)
}

// 添加新的公共函数用于单个文件的查询（保持原有接口不变）
pub fn query_single_file(sample_path: &str, db_path: &str, min_ani: f64, mismatch: usize, min_tags_genome: usize) -> Result<Vec<QueryResult>> {
    eprintln!("Processing database file: {}", db_path);
    
    // 读取数据库文件
    let db_file = File::open(db_path)
        .with_context(|| format!("Failed to open database file: {}", db_path))?;
    let db_reader = BufReader::new(db_file);
    let db_entries: Vec<SyldbEntry> = bincode::deserialize_from(db_reader)
        .with_context(|| format!("Failed to deserialize database file: {}", db_path))?;

    eprintln!("Found {} entries in database", db_entries.len());
    if let Some(enzyme) = detect_db_enzyme(&db_entries) {
        eprintln!("Database enzyme: {}", enzyme);
    }

    // 读取样本文件 - 优化大文件读取
    let sample_file = File::open(sample_path)
        .with_context(|| format!("Failed to open sample file: {}", sample_path))?;
    let sample_reader = BufReader::with_capacity(100_000_000, sample_file); // 100MB 缓冲区
    let sample_entries: Vec<SylspEntry> = bincode::deserialize_from(sample_reader)
        .with_context(|| format!("Failed to deserialize sample file: {}", sample_path))?;

    eprintln!("Found {} entries in sample", sample_entries.len());

    // 检查样本数据的有效性
    if sample_entries.is_empty() {
        eprintln!("Warning: Sample {} has no tags", sample_path);
        return Ok(Vec::new());
    }

    // 酶集合一致性检查（无 CLI 上下文：严格模式）
    check_enzyme_compat(detect_db_enzyme(&db_entries).as_deref(), detect_sample_enzyme(&sample_entries).as_deref(), false, sample_path)?;

    // 按样本源分组
    let mut sample_groups: FxHashMap<String, Vec<&SylspEntry>> = FxHashMap::default();
    for entry in &sample_entries {
        sample_groups.entry(entry.sample_source.clone())
            .or_default()
            .push(entry);
    }

    // 并行处理每个样本组，然后合并结果
    let all_results: Vec<QueryResult> = sample_groups.par_iter()
        .flat_map(|(sample_source, entries)| {
            // 统计样本中每个 tag 的覆盖度（重数）
            let sample_counts = sample_tag_counts_ref(entries);
            let total_sample_tags = entries.len();
            // 无 CLI 上下文：mm1 默认使用数据驱动的 e 估计
            let error_rate = sample_error_rate(&db_entries, &sample_counts, mismatch, None, false);

            // 并行处理每个基因组进行比对
            db_entries.par_iter().filter_map(|db_entry| {
                let covs: Vec<u32> = db_entry.tags.iter().enumerate()
                    .filter_map(|(i, t)| {
                        let seq = db_entry.tag_seqs.as_ref().and_then(|v| v.get(i));
                        lookup_tag_coverage(*t, seq, &sample_counts, mismatch)
                    })
                    .collect();

                let mut result = calculate_statistics(covs, &db_entry.tag_lengths, total_sample_tags, db_entry.tags.len(), mismatch, error_rate);
                result.sample_file = sample_source.clone();
                result.genome_file = db_path.to_string();
                result.contig_name = db_entry.sequence_id.clone();

                if filter_results(&result, Some(min_ani), min_tags_genome) {
                    Some(result)
                } else {
                    None
                }
            }).collect::<Vec<QueryResult>>()
        })
        .collect();

    Ok(all_results)
}

// 从文件路径中提取基因组ID的函数
// 基因组 ID 与输入文件名逐字一致（不剥离扩展名）
fn extract_genome_id(file_path: &str) -> String {
    let path = std::path::Path::new(file_path);
    
    // 获取文件名（不含路径）
    let file_name = path.file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("");
    
    file_name.to_string()
}

// 读取taxonomy文件并建立genome到分类信息的映射
fn read_taxonomy_file(taxonomy_file: &str) -> Result<FxHashMap<String, Arc<TaxonomyInfo>>> {
    use std::io::BufRead;
    
    let file = File::open(taxonomy_file)
        .with_context(|| format!("Failed to open taxonomy file: {}", taxonomy_file))?;
    let reader = BufReader::new(file);
    
    // 预分配 HashMap 容量以提高性能
    let mut taxonomy_map = FxHashMap::default(); // 使用默认容量
    let mut line_count = 0;
    
    for line in reader.lines() {
        line_count += 1;
        let line = line.context("Failed to read line from taxonomy file")?;
        let line = line.trim();
        
        // 跳过标题行和空行
        if line_count == 1 || line.is_empty() || line.starts_with('#') {
            continue;
        }
        
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 3 {
            eprintln!("Warning: Invalid line format at line {}: {}", line_count, line);
            continue;
        }
        
        let accession = parts[0].trim();
        let gtdb_taxonomy = parts[1].trim();
        
        // 处理genome ID：去掉前缀 (RS_, GB_, 等)
        let genome_id = if accession.starts_with("RS_") {
            &accession[3..]
        } else if accession.starts_with("GB_") {
            &accession[3..]
        } else if accession.starts_with("GS_") {
            &accession[3..]
        } else {
            accession
        };
        
        // 解析GTDB分类信息并包装在 Arc 中
        match TaxonomyInfo::from_gtdb_string(gtdb_taxonomy) {
            Ok(taxonomy) => {
                let taxonomy_arc = Arc::new(taxonomy);
                // 添加基本ID（例如：GCF_000006685.1）
                taxonomy_map.insert(genome_id.to_string(), Arc::clone(&taxonomy_arc));
                
                // 添加带_genomic后缀的ID（例如：GCF_000006685.1_genomic）
                let genomic_id = format!("{}_genomic", genome_id);
                taxonomy_map.insert(genomic_id, taxonomy_arc);
            }
            Err(e) => {
                eprintln!("Warning: Failed to parse taxonomy for {}: {}", genome_id, e);
            }
        }
    }
    
    eprintln!("Loaded taxonomy information for {} genome variants", taxonomy_map.len());
    Ok(taxonomy_map)
}

// 从基因组级别聚合到物种级别
fn aggregate_to_species_level(
    sample_groups: &HashMap<String, Vec<GenomeProfileResult>>,
    taxonomy_map: &FxHashMap<String, Arc<TaxonomyInfo>>,
    effective_min_ani: f64,
) -> Result<Vec<SpeciesAbundanceResult>> {
    use std::sync::Mutex;
    
    let species_map = Arc::new(Mutex::new(FxHashMap::<String, SpeciesAbundanceResult>::default()));
    
    // 采用 sylph 的分层并行策略进行物种聚合
    let sample_groups_arc = Arc::new(sample_groups);
    
    // 外层并行：处理样本组
    sample_groups_arc.par_iter()
        .for_each(|(sample_id, genome_results)| {
            // 为每个样本组创建局部聚合结果 - 预分配容量
            let mut local_species_map: FxHashMap<String, SpeciesAbundanceResult> = FxHashMap::default();
            
            for genome_result in genome_results {
                // 额外的过滤条件：确保只有高质量的genome参与物种聚合
                if genome_result.adjusted_ani < effective_min_ani ||
                   genome_result.common_tags == 0 {
                    continue;
                }
                
                // 从genome_id中提取标准化的genome标识符
                let genome_id = extract_genome_id_from_path(&genome_result.genome_id);
                
                // 查找对应的分类信息：先按逐字 ID 查，再回退到剥离扩展名的旧式 ID
                // （taxonomy 文件的键历史上不带 .fasta 等扩展名）
                let taxonomy_arc = taxonomy_map.get(genome_id).or_else(|| {
                    strip_genome_extension(genome_id).and_then(|stripped| taxonomy_map.get(stripped))
                });
                if let Some(taxonomy_arc) = taxonomy_arc {
                    let species_key = taxonomy_arc.get_species_key();
                    
                    // 获取或创建物种条目 - 使用 Arc 共享而不是克隆
                    let species_result = local_species_map.entry(species_key).or_insert_with(|| {
                        SpeciesAbundanceResult {
                            taxonomy: Arc::clone(taxonomy_arc),
                            sample_abundances: FxHashMap::default(),
                            total_tags: 0,
                            genome_count: 0,
                            reads_count: 0,
                            gscore: 0.0,
                        }
                    });
                    
                    // 累加样本丰度
                    *species_result.sample_abundances.entry(sample_id.clone()).or_insert(0.0) += 
                        genome_result.taxonomic_abundance;
                    
                    // 累加标签数、基因组计数和 reads 数
                    species_result.total_tags += genome_result.total_tags;
                    species_result.genome_count += 1;
                    // 使用 common_tags 作为该基因组在该样本中的 reads 数代理
                    species_result.reads_count += genome_result.common_tags;
                } else {
                    eprintln!("Warning: No taxonomy information found for genome: {}", genome_id);
                }
            }
            
            // 将局部结果合并到全局结果中
            let mut global_map = species_map.lock().unwrap();
            for (species_key, local_result) in local_species_map {
                let global_result = global_map.entry(species_key).or_insert_with(|| {
                    SpeciesAbundanceResult {
                        taxonomy: Arc::clone(&local_result.taxonomy),
                        sample_abundances: FxHashMap::default(),
                        total_tags: 0,
                        genome_count: 0,
                        reads_count: 0,
                        gscore: 0.0,
                    }
                });
                
                // 合并样本丰度
                for (sample_id, abundance) in local_result.sample_abundances {
                    *global_result.sample_abundances.entry(sample_id).or_insert(0.0) += abundance;
                }
                
                // 合并标签数、基因组计数和 reads 数
                global_result.total_tags += local_result.total_tags;
                global_result.genome_count += local_result.genome_count;
                global_result.reads_count += local_result.reads_count;
            }
        });
    
    let species_map = Arc::try_unwrap(species_map).unwrap().into_inner().unwrap();
    let mut results: Vec<SpeciesAbundanceResult> = species_map.into_values().collect();
    
    // 按分类学层级排序
    results.sort_by(|a, b| {
        a.taxonomy.kingdom.cmp(&b.taxonomy.kingdom)
            .then_with(|| a.taxonomy.phylum.cmp(&b.taxonomy.phylum))
            .then_with(|| a.taxonomy.class.cmp(&b.taxonomy.class))
            .then_with(|| a.taxonomy.order.cmp(&b.taxonomy.order))
            .then_with(|| a.taxonomy.family.cmp(&b.taxonomy.family))
            .then_with(|| a.taxonomy.genus.cmp(&b.taxonomy.genus))
            .then_with(|| a.taxonomy.species.cmp(&b.taxonomy.species))
    });
    
    Ok(results)
}

// 从文件路径或genome_id中提取标准化的genome标识符
// 注意：不剥离扩展名，基因组 ID 与输入文件名逐字一致（与 sylph 行为相同），
// 保证下游 ground truth / 评分脚本能直接按文件名 join。
fn extract_genome_id_from_path(input: &str) -> &str {
    // 如果输入包含路径分隔符，提取文件名
    if input.contains('/') {
        std::path::Path::new(input)
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or(input)
    } else {
        input
    }
}

// 旧式 ID 形式：剥离已知基因组文件扩展名（.fasta/.fa/.fna 及其 .gz 变体）。
// 仅用于 taxonomy 查询的兼容回退（taxonomy 文件的键历史上不带扩展名）。
fn strip_genome_extension(id: &str) -> Option<&str> {
    id.strip_suffix(".fasta.gz").or_else(|| id.strip_suffix(".fasta"))
        .or_else(|| id.strip_suffix(".fa.gz")).or_else(|| id.strip_suffix(".fa"))
        .or_else(|| id.strip_suffix(".fna.gz")).or_else(|| id.strip_suffix(".fna"))
}

// 从db文件中读取基因组映射关系
fn read_genome_mapping(db_path: &str) -> Result<FxHashMap<String, (String, String)>> {
    let db_file = File::open(db_path)?;
    let db_reader = BufReader::new(db_file);
    let db_entries: Vec<SyldbEntry> = bincode::deserialize_from(db_reader)
        .with_context(|| format!("Failed to deserialize database file: {}", db_path))?;
    
    // 并行处理数据库条目生成映射
    let genome_map: FxHashMap<String, (String, String)> = db_entries.par_iter()
        .map(|entry| {
            // 获取原始基因组文件路径
            let genome_source = entry.genome_source.clone();
            let genome_id = if let Some(file_name) = std::path::Path::new(&genome_source)
                .file_name()
                .and_then(|s| s.to_str()) 
            {
                // 基因组 ID 与输入文件名逐字一致（不剥离扩展名）
                file_name.to_string()
            } else {
                genome_source.clone()
            };
            
            // 返回(序列ID, (基因组ID, 基因组源文件))
            (entry.sequence_id.clone(), (genome_id, genome_source))
        })
        .collect();
    
    Ok(genome_map)
}

// 读取样本文件列表
fn read_sample_list(list_file: &str) -> Result<Vec<String>> {
    let content = std::fs::read_to_string(list_file)?;
    // 预分配 Vec 容量并减少字符串克隆
    let lines: Vec<String> = content.lines()
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
        .map(|s| s.to_string())
        .collect();
    Ok(lines)
}

// 生成TSV格式的丰度矩阵
fn write_abundance_matrix(
    sample_groups: &HashMap<String, Vec<GenomeProfileResult>>,
    all_genomes: &HashSet<String>,
    log_path: Option<String>,
    tsv_name: &str,
    writer: &mut Box<dyn Write + Send>,
) -> Result<()> {
    // 如果指定了log_path，使用它，否则使用当前目录
    let output_dir = if let Some(path) = log_path {
        PathBuf::from(path)
    } else {
        std::env::current_dir()?
    };

    // 确保目录存在
    std::fs::create_dir_all(&output_dir)?;

    // 构建TSV文件路径
    let tsv_path = output_dir.join(tsv_name);
    let mut tsv_writer = BufWriter::new(File::create(tsv_path)?);

    // 获取所有样本ID并排序
    let mut sample_ids: Vec<_> = sample_groups.keys().collect();
    sample_ids.sort();

    // 写入表头
    write!(tsv_writer, "Genome")?;
    write!(writer, "\nAbundance Matrix:\n")?;
    write!(writer, "Genome")?;
    for sample_id in &sample_ids {
        write!(tsv_writer, "\t{}", sample_id)?;
        write!(writer, "\t{}", sample_id)?;
    }
    writeln!(tsv_writer)?;
    writeln!(writer)?;

    // 采用 sylph 的高效并行数据收集策略
    let genome_data: Vec<(String, Vec<f64>)> = all_genomes.par_iter()
        .map(|genome_id| {
            let abundances: Vec<f64> = sample_ids.iter()
                .map(|sample_id| {
                    sample_groups.get(sample_id.as_str())
                        .and_then(|results| results.iter()
                            .find(|r| r.genome_id == *genome_id))
                        .map(|r| r.taxonomic_abundance)
                        .unwrap_or(0.0)
                })
                .collect();
            (genome_id.clone(), abundances)
        })
        .collect();

    // 写入每个基因组的丰度数据
    for (genome_id, abundances) in genome_data {
        write!(tsv_writer, "{}", genome_id)?;
        write!(writer, "{}", genome_id)?;
        for abundance in abundances {
            write!(tsv_writer, "\t{:.4}", abundance)?;
            write!(writer, "\t{:.4}", abundance)?;
        }
        writeln!(tsv_writer)?;
        writeln!(writer)?;
    }
    writeln!(writer)?;

    Ok(())
}

// 生成物种级别的TSV格式丰度矩阵
fn write_species_abundance_matrix(
    species_results: &[SpeciesAbundanceResult],
    all_samples: &HashSet<String>,
    log_path: Option<String>,
    tsv_name: &str,
    writer: &mut Box<dyn Write + Send>,
) -> Result<()> {
    // 如果指定了log_path，使用它，否则使用当前目录
    let output_dir = if let Some(path) = log_path {
        PathBuf::from(path)
    } else {
        std::env::current_dir()?
    };

    // 确保目录存在
    std::fs::create_dir_all(&output_dir)?;

    // 构建TSV文件路径
    let tsv_path = output_dir.join(tsv_name);
    let mut tsv_writer = BufWriter::new(File::create(tsv_path)?);

    // 获取所有样本ID并排序
    let mut sample_ids: Vec<_> = all_samples.iter().collect();
    sample_ids.sort();

    // 写入表头 (参考Abundance_Stat.all.xls格式)
    write!(tsv_writer, "#Kingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies")?;
    write!(writer, "\nSpecies-level Abundance Matrix:\n")?;
    write!(writer, "#Kingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies")?;
    for sample_id in &sample_ids {
        write!(tsv_writer, "\t{}", sample_id)?;
        write!(writer, "\t{}", sample_id)?;
    }
    writeln!(tsv_writer)?;
    writeln!(writer)?;

    // 采用 sylph 的高效并行数据收集策略
    let species_data: Vec<(Arc<TaxonomyInfo>, Vec<f64>)> = species_results.par_iter()
        .map(|species_result| {
            let abundances: Vec<f64> = sample_ids.iter()
                .map(|sample_id| {
                    species_result.sample_abundances
                        .get(sample_id.as_str())
                        .copied()
                        .unwrap_or(0.0)
                })
                .collect();
            (Arc::clone(&species_result.taxonomy), abundances)
        })
        .collect();

    // 写入每个物种的丰度数据
    for (taxonomy_arc, abundances) in species_data {
        // 写入分类学信息（7列）
        write!(tsv_writer, "{}\t{}\t{}\t{}\t{}\t{}\t{}", 
               taxonomy_arc.kingdom, taxonomy_arc.phylum, taxonomy_arc.class,
               taxonomy_arc.order, taxonomy_arc.family, taxonomy_arc.genus, taxonomy_arc.species)?;
        write!(writer, "{}\t{}\t{}\t{}\t{}\t{}\t{}", 
               taxonomy_arc.kingdom, taxonomy_arc.phylum, taxonomy_arc.class,
               taxonomy_arc.order, taxonomy_arc.family, taxonomy_arc.genus, taxonomy_arc.species)?;
        
        // 写入各个样本的丰度值
        for abundance in abundances {
            write!(tsv_writer, "\t{:.6}", abundance)?;
            write!(writer, "\t{:.6}", abundance)?;
        }
        writeln!(tsv_writer)?;
        writeln!(writer)?;
    }
    writeln!(writer)?;

    Ok(())
}

// 从缓存的数据库条目中构建基因组映射关系
fn build_genome_mapping_from_cache(cached_db_entries: &[SyldbEntry]) -> FxHashMap<String, (String, String)> {
    // 预分配 HashMap 容量以提高性能
    let mut genome_map = FxHashMap::default();
    
    // 并行处理数据库条目生成映射
    for entry in cached_db_entries {
        // 获取原始基因组文件路径
        let genome_source = entry.genome_source.clone();
        let genome_id = if let Some(file_name) = std::path::Path::new(&genome_source)
            .file_name()
            .and_then(|s| s.to_str()) 
        {
            // 基因组 ID 与输入文件名逐字一致（不剥离扩展名）
            file_name.to_string()
        } else {
            genome_source.clone()
        };
        
        // 返回(序列ID, (基因组ID, 基因组源文件))
        genome_map.insert(entry.sequence_id.clone(), (genome_id, genome_source));
    }
    
    genome_map
}

// 更新profile函数
pub fn profile(args: ProfileArgs) -> Result<()> {
    // sketch 模式下 mismatch 需要 tag 序列生成 1-mismatch 邻居，而 k-mer sketch 不存序列。
    // 该检查必须在任何文件 IO / 线程池初始化之前，保证快速失败。
    if args.sketch_mode && args.mismatch != 0 {
        return Err(anyhow!(
            "--sketch-mode requires --mismatch 0: k-mer sketches (.db/.sp) do not store tag sequences, so 1-mismatch neighbor generation is impossible"
        ));
    }

    // 处理minimum_ani参数：如果没有传入参数，使用默认值
    let effective_min_ani = args.minimum_ani.unwrap_or(PROFILE_MIN_ANI);
    let min_eff_coverage = args.min_eff_coverage;
    let min_shared_tags = args.min_shared_tags;
    let min_tags_genome = args.min_tags_genome;
    eprintln!("Using minimum ANI threshold: {:.1}%", effective_min_ani);
    eprintln!("Using minimum shared tags threshold: {}", min_shared_tags);
    eprintln!("Using minimum tags per genome threshold: {}", min_tags_genome);
    if min_eff_coverage > 0.0 {
        eprintln!("Using minimum effective coverage threshold: {:.2}", min_eff_coverage);
    }
    
    // 优化线程池配置 - 采用 sylph 的策略
    let _max_ram = args.threads * 2; // 简单的内存限制，每线程2GB
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .context("Failed to initialize thread pool")?;

    // 一次性读取并缓存数据库文件 - 优化大文件读取
    eprintln!("Loading database file: {}", args.db_file);
    
    let db_file = File::open(&args.db_file)
        .with_context(|| format!("Failed to open database file: {}", args.db_file))?;
    let db_reader = BufReader::with_capacity(100_000_000, db_file); // 100MB 缓冲区
    let cached_db_entries: Vec<SyldbEntry> = if args.sketch_mode {
        // sketch 模式：.db 是 Vec<GenomeSketch>，逐条适配成 SyldbEntry
        let genome_sketches: Vec<GenomeSketch> = bincode::deserialize_from(db_reader)
            .with_context(|| format!("Failed to deserialize sketch database file (--sketch-mode expects Vec<GenomeSketch>): {}", args.db_file))?;
        genome_sketches
            .iter()
            .map(genome_sketch_to_syldb_entry)
            .collect()
    } else {
        bincode::deserialize_from(db_reader)
            .with_context(|| format!("Failed to deserialize database file: {}", args.db_file))?
    };
    // 关键：按基因组聚合所有 contig（去重 tag），整个 profile 流程按"整基因组"统计
    let cached_db_entries = aggregate_db_by_genome(cached_db_entries);

    eprintln!("Cached {} genomes from database (after aggregating contigs)", cached_db_entries.len());
    let db_enzyme = detect_db_enzyme(&cached_db_entries);
    if let Some(enzyme) = &db_enzyme {
        eprintln!("Database enzyme: {}", enzyme);
    }
    ensure_tag_seqs_for_mismatch(&cached_db_entries, args.mismatch)?;

    // 一次性读取并缓存所有样本文件 - 优化大文件读取
    eprintln!("Loading sample files: {}", args.sample_file);
    let sample_files: Vec<String> = if args.sample_file.ends_with(".txt") {
        read_sample_list(&args.sample_file)?
    } else {
        vec![args.sample_file.clone()]
    };

    let mut cached_sample_entries: FxHashMap<String, Vec<SylspEntry>> = FxHashMap::default();
    let mut cached_sample_sketches: FxHashMap<String, Vec<SequencesSketch>> = FxHashMap::default();
    for sample_path in &sample_files {
        if args.sketch_mode {
            // sketch 模式：.sp 是 SequencesSketch（或合并的 Vec<SequencesSketch>）。
            // kmer_counts 保持 map 形式驻留，绝不展开成逐条观测（大样本会爆内存）。
            let sketches = load_sequences_sketches(sample_path)?;
            let total_kmers: usize = sketches.iter().map(|s| s.kmer_counts.len()).sum();
            eprintln!(
                "Cached {} sketch(es) with {} distinct k-mers from {}",
                sketches.len(),
                total_kmers,
                sample_path
            );
            cached_sample_sketches.insert(sample_path.clone(), sketches);
        } else {
            let sample_file = File::open(sample_path)
                .with_context(|| format!("Failed to open sample file: {}", sample_path))?;
            let sample_reader = BufReader::with_capacity(100_000_000, sample_file); // 100MB 缓冲区
            let sample_entries: Vec<SylspEntry> = bincode::deserialize_from(sample_reader)
                .with_context(|| format!("Failed to deserialize sample file: {}", sample_path))?;
            cached_sample_entries.insert(sample_path.clone(), sample_entries);
        }
    }
    eprintln!("Cached {} sample files", sample_files.len());

    // 酶集合一致性检查：在并行流程开始前，对每个样本 sketch 与 DB 比较
    // （query_single_file_with_cached_db 的 Result 会被 `if let Ok` 吞掉，硬错误必须放在这里）。
    // sketch 模式下两侧都用 sketch_enzyme_tag 生成自描述串，c/k 一致即天然通过。
    for sample_path in &sample_files {
        if args.sketch_mode {
            if let Some(sketches) = cached_sample_sketches.get(sample_path) {
                for sketch in sketches {
                    let sample_enzyme = sketch_enzyme_tag(sketch.c, sketch.k);
                    check_enzyme_compat(
                        db_enzyme.as_deref(),
                        Some(&sample_enzyme),
                        args.allow_enzyme_mismatch,
                        sample_path,
                    )?;
                }
            }
        } else if let Some(entries) = cached_sample_entries.get(sample_path) {
            check_enzyme_compat(
                db_enzyme.as_deref(),
                detect_sample_enzyme(entries).as_deref(),
                args.allow_enzyme_mismatch,
                sample_path,
            )?;
        }
    }

    // 从缓存的数据库构建基因组映射关系
    let genome_mapping = build_genome_mapping_from_cache(&cached_db_entries);

    // 创建输出写入器
    let mut writer = create_multi_writer(&args.out_file_name)?;

    // 存储所有样本的结果 - 预分配容量，使用 Mutex 保护
    let all_results = Arc::new(Mutex::new(
        FxHashMap::<(String, String), GenomeProfileResult>::default(),
    ));

    // 采用 sylph 的简化并行处理策略
    let step = usize::max(
        args.threads / 3 + 1,
        usize::min(sample_files.len(), args.threads),
    );
    let chunks: Vec<Vec<String>> = sample_files
        .chunks(step)
        .map(|chunk| chunk.to_vec())
        .collect();

    // 使用 sylph 风格的分块处理，集成k-mer重新分配机制
    chunks.into_iter().for_each(|chunk| {
        chunk.into_par_iter().for_each(|sample_file| {
            // 第一阶段：计算初步结果（不使用重新分配）
            let initial_results_opt: Option<Vec<QueryResult>> = if args.sketch_mode {
                // sketch 模式：每个 SequencesSketch 是一个样本，直接借用 kmer_counts，
                // 不展开成逐条观测。
                cached_sample_sketches.get(&sample_file).map(|sketches| {
                    let results: Vec<QueryResult> = sketches
                        .par_iter()
                        .filter(|s| !s.kmer_counts.is_empty())
                        .flat_map(|sketch| {
                            let (sample_source, sample_counts, total_sample_tags) =
                                sketch_sample_view(sketch);
                            eprintln!(
                                "Processing sample sketch: {} with {} distinct k-mers",
                                sample_source,
                                sample_counts.len()
                            );
                            query_one_sample_with_cached_db(
                                &sample_source,
                                sample_counts,
                                total_sample_tags,
                                &args.db_file,
                                &cached_db_entries,
                                effective_min_ani,
                                args.mismatch,
                                min_shared_tags,
                                min_tags_genome,
                                args.read_error_rate,
                                args.no_error_correction,
                            )
                        })
                        .collect();
                    results
                })
            } else {
                query_single_file_with_cached_db(
                    &sample_file,
                    &args.db_file,
                    &cached_db_entries,
                    &cached_sample_entries,
                    effective_min_ani,
                    args.mismatch,
                    min_shared_tags,
                    min_tags_genome,
                    args.read_error_rate,
                    args.no_error_correction,
                )
                .ok()
            };

            if let Some(initial_results) = initial_results_opt {
                // 按ANI排序
                let mut initial_results = initial_results;
                initial_results
                    .sort_by(|a, b| b.adjusted_ani.partial_cmp(&a.adjusted_ani).unwrap());

                // 第二阶段：构建winner table并重新分配（模仿sylph的两阶段处理）
                eprintln!(
                    "{} taxonomic profiling; reassigning tags for {} genomes...",
                    &sample_file,
                    initial_results.len()
                );

                // 构建winner table
                let (winner_map, interner) =
                    build_winner_table(&initial_results, &cached_db_entries, true); // 启用日志

                // 使用winner table重新计算结果。
                // sketch 模式下没有 SylspEntry 列表；recalculate_abundances_after_reassignment
                // 只用结果、不读样本条目，传空切片即可。
                let empty_sample_entries: &[SylspEntry] = &[];
                let sample_entries_opt: Option<&[SylspEntry]> = if args.sketch_mode {
                    Some(empty_sample_entries)
                } else {
                    cached_sample_entries
                        .get(&sample_file)
                        .map(|v| v.as_slice())
                };
                if let Some(sample_entries) = sample_entries_opt {
                    let mut reassigned_results = if args.sketch_mode {
                        let sketches = cached_sample_sketches
                            .get(&sample_file)
                            .expect("sample sketch cache populated during loading");
                        // 初始结果索引（与 recalculate_with_winner_table 内部一致）
                        let initial_by_contig: FxHashMap<&str, &QueryResult> = initial_results
                            .iter()
                            .map(|r| (r.contig_name.as_str(), r))
                            .collect();
                        let mut out = Vec::new();
                        for sketch in sketches {
                            if sketch.kmer_counts.is_empty() {
                                continue;
                            }
                            let (sample_source, sample_counts, total_sample_tags) =
                                sketch_sample_view(sketch);
                            let mut res = recalculate_one_sample_with_winner_table(
                                &sample_source,
                                sample_counts,
                                total_sample_tags,
                                &cached_db_entries,
                                &winner_map,
                                &interner,
                                effective_min_ani,
                                false,
                                args.mismatch,
                                min_shared_tags,
                                min_tags_genome,
                                args.read_error_rate,
                                args.no_error_correction,
                                &initial_by_contig,
                                args.reassign_protection,
                            );
                            out.append(&mut res);
                        }
                        out
                    } else {
                        recalculate_with_winner_table(
                            &cached_db_entries,
                            sample_entries,
                            &winner_map,
                            &interner,
                            effective_min_ani,
                            false,
                            args.mismatch,
                            min_shared_tags,
                            min_tags_genome,
                            args.read_error_rate,
                            args.no_error_correction,
                            &initial_results,
                            args.reassign_protection,
                        )
                    };

                    // 第三阶段：过滤过度重新分配的基因组
                    reassigned_results = filter_over_reassigned_genomes(
                        &initial_results,
                        &reassigned_results,
                        effective_min_ani,
                        K,
                    );

                    // 第四阶段：重新计算丰度
                    recalculate_abundances_after_reassignment(
                        &mut reassigned_results,
                        sample_entries,
                    );

                    eprintln!(
                        "{} has {} genomes passing profiling threshold after reassignment.",
                        &sample_file,
                        reassigned_results.len()
                    );

                    // 按基因组ID分组结果 - 修复：确保每个样本源都被正确处理
                    for result in reassigned_results {
                        if let Some((genome_id, _)) = genome_mapping.get(&result.contig_name) {
                            // 关键修复：使用实际的样本源ID作为key的一部分
                            let key = (genome_id.clone(), result.sample_file.clone());
                            let mut all_results = all_results.lock().unwrap();
                            let entry = all_results.entry(key).or_insert_with(|| {
                                GenomeProfileResult {
                                    genome_id: genome_id.clone(),
                                    sample_id: result.sample_file.clone(), // 这里保存的是实际的样本源ID
                                    file_path: sample_file.clone(),
                                    adjusted_ani: 0.0,
                                    taxonomic_abundance: 0.0,
                                    sequence_abundance: 0.0,
                                    common_tags: 0,
                                    total_tags: 0,
                                    eff_cov: 0.0,
                                }
                            });

                            // 累加标签数
                            entry.common_tags += result.shared_tags;
                            entry.total_tags += result.ref_tags;
                            entry.eff_cov += result.eff_cov;

                            // 使用共享标签数作为权重计算加权平均ANI
                            if entry.common_tags > 0 {
                                entry.adjusted_ani = (entry.adjusted_ani
                                    * (entry.common_tags - result.shared_tags) as f64
                                    + result.adjusted_ani * result.shared_tags as f64)
                                    / entry.common_tags as f64;
                            }
                        }
                    }
                }
            }
        });
    });
    
    // 收集所有基因组ID
    let mut all_genomes: HashSet<String> = HashSet::new();
    for entry in genome_mapping.values() {
        all_genomes.insert(entry.0.clone());
    }

    // 转换为向量以便排序和分组
    let results: Vec<_> = all_results.lock().unwrap().values().cloned().collect();
    
    // 按样本分组计算丰度
    let mut sample_groups: HashMap<String, Vec<GenomeProfileResult>> = HashMap::new();
    for result in results {
        // 使用实际的样本来源而不是文件名
        sample_groups.entry(result.sample_id.clone())
            .or_default()
            .push(result);
    }
    
    // 采用 sylph 的简单策略 - 顺序计算丰度，避免复杂的并行迭代器组合
    for (_sample_id, group) in sample_groups.iter_mut() {
        // 按ANI排序（参考sylph的排序机制）
        group.sort_by(|a, b| b.adjusted_ani.partial_cmp(&a.adjusted_ani).unwrap());
        
        // 过滤掉不符合profile要求的genome（依赖覆盖度校正后的 ANI + 基因组大小 + 最小共享 tag 数 + effective coverage）
        group.retain(|r| {
            r.common_tags > 0 &&
            r.common_tags >= min_shared_tags &&
            r.adjusted_ani >= effective_min_ani &&
            r.total_tags >= min_tags_genome &&
            r.eff_cov >= min_eff_coverage
        });
        
        // 计算总覆盖度，包括所有检测到的标签
        let total_genome_cov: f64 = group.iter()
            .map(|r| if r.common_tags > 0 { r.eff_cov } else { 0.0 })
            .sum();
        
        let total_seq_cov: f64 = group.iter()
            .map(|r| if r.common_tags > 0 { 
                r.eff_cov * r.total_tags as f64 
            } else { 
                0.0 
            })
            .sum();
        
        // 计算每个结果的丰度 - 采用 sylph 的顺序处理方式
        for result in group.iter_mut() {
            // 只要有共享标签就计算丰度
            if result.common_tags > 0 {
                result.taxonomic_abundance = if total_genome_cov > 0.0 {
                    result.eff_cov / total_genome_cov * 100.0
                } else {
                    0.0
                };
                
                result.sequence_abundance = if total_seq_cov > 0.0 {
                    result.eff_cov * result.total_tags as f64 / total_seq_cov * 100.0
                } else {
                    0.0
                };
            } else {
                result.taxonomic_abundance = 0.0;
                result.sequence_abundance = 0.0;
            }
        }
    }

    // 检查是否提供了taxonomy文件以进行物种级别聚合
    if let Some(taxonomy_file) = &args.taxonomy_file {
        eprintln!("Loading taxonomy information from: {}", taxonomy_file);
        
        // 读取分类学信息
        let taxonomy_map = read_taxonomy_file(taxonomy_file)?;
        
        // 聚合到物种级别
        let mut species_results = aggregate_to_species_level(&sample_groups, &taxonomy_map, effective_min_ani)?;
        
        // 获取所有样本ID
        let all_samples: HashSet<String> = sample_groups.keys().cloned().collect();
        
        // 生成过滤前的物种级别TSV格式丰度矩阵
        let pre_filter_tsv_name = format!("pre_gscore_filter_{}", args.tsv_name);
        eprintln!("Writing pre-filter species abundance matrix: {}", pre_filter_tsv_name);
        write_species_abundance_matrix(&species_results, &all_samples, args.log_path.clone(), &pre_filter_tsv_name, &mut writer)?;
        
        // 应用 G-score 过滤
        eprintln!("Applying G-score filtering with threshold: {:.2}", args.gscore_threshold);
        species_results = filter_species_by_gscore(&mut species_results, args.gscore_threshold);
        
        // 生成过滤后的物种级别TSV格式丰度矩阵
        eprintln!("Writing post-filter species abundance matrix: {}", args.tsv_name);
        write_species_abundance_matrix(&species_results, &all_samples, args.log_path.clone(), &args.tsv_name, &mut writer)?;
        
        // 输出物种级别的统计信息
        writeln!(writer, "Species-level Profile Results:")?;
        writeln!(writer, "------------------------------")?;
        writeln!(writer, "Sample files: {} files processed", sample_files.len())?;
        writeln!(writer, "Database file: {}", args.db_file)?;
        writeln!(writer, "Taxonomy file: {}", taxonomy_file)?;
        writeln!(writer, "Total species detected: {}", species_results.len())?;
        writeln!(writer, "\nSpecies composition summary:")?;
        writeln!(writer, "{:<50} {:<15} {:<15} {:<15} {:<10} {:<10}", 
            "Species", "Genomes", "Total_Tags", "Reads_Count", "G-score", "Avg_Abundance")?;
        writeln!(writer, "{:-<120}", "")?;
        
        for species_result in &species_results {
            let avg_abundance: f64 = species_result.sample_abundances.values().sum::<f64>() 
                / species_result.sample_abundances.len() as f64;
            let species_name = if species_result.taxonomy.species.is_empty() {
                format!("{}_sp", species_result.taxonomy.genus)
            } else {
                species_result.taxonomy.species.clone()
            };
            
            writeln!(writer, "{:<50} {:<15} {:<15} {:<15} {:<10.2} {:<10.4}", 
                species_name,
                species_result.genome_count,
                species_result.total_tags,
                species_result.reads_count,
                species_result.gscore,
                avg_abundance)?;
        }
        
    } else {
        // 原始的基因组级别输出
        // 生成TSV格式的丰度矩阵
        write_abundance_matrix(&sample_groups, &all_genomes, args.log_path.clone(), &args.tsv_name, &mut writer)?;

        // 将所有结果收集到一个新的向量中
        let mut final_results: Vec<GenomeProfileResult> = sample_groups.into_values().flatten().collect();
        
        // 按基因组ID和ANI值排序
        final_results.sort_by(|a, b| {
            a.genome_id.cmp(&b.genome_id)
                .then_with(|| b.adjusted_ani.partial_cmp(&a.adjusted_ani).unwrap())
        });
        
        // 输出结果
        writeln!(writer, "Genome-level Profile Results:")?;
        writeln!(writer, "-----------------------------")?;
        writeln!(writer, "Sample files: {} files processed", sample_files.len())?;
        writeln!(writer, "Database file: {}", args.db_file)?;
        writeln!(writer, "\nGenome composition:")?;
        writeln!(writer, "{:<30} {:<20} {:<10} {:<12} {:<12} {:<12} {:<12} {:<10}", 
            "Genome_ID", "Sample_ID", "ANI(%)", "Tax_Abund(%)", "Seq_Abund(%)", "Common_Tags", "Total_Tags", "Eff_cov")?;
        writeln!(writer, "{:-<110}", "")?;
        
        let mut current_genome = String::new();
        for result in final_results {
            if current_genome != result.genome_id {
                if !current_genome.is_empty() {
                    writeln!(writer)?;
                }
                current_genome = result.genome_id.clone();
            }
            
            writeln!(writer, "{:<30} {:<20} {:<10.2} {:<12.2} {:<12.2} {:<12} {:<12} {:<10.3}", 
                result.genome_id,
                result.sample_id,  // 使用实际的样本来源
                result.adjusted_ani,
                result.taxonomic_abundance,
                result.sequence_abundance,
                result.common_tags,
                result.total_tags,
                result.eff_cov)?;
        }
    }
    
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_genome_sketch_adapter_fields() {
        // Adapter must carry the k-mers over as tags, set every tag length to k,
        // derive a self-describing enzyme string, and leave sequence-bearing fields empty.
        let sketch = GenomeSketch {
            file_name: "genome1".to_string(),
            first_contig_name: "contig1".to_string(),
            gn_size: 100_000,
            c: 200,
            k: 31,
            min_spacing: 0,
            genome_kmers: vec![5, 1, 9, 1, 5, 42],
        };
        let entry = genome_sketch_to_syldb_entry(&sketch);
        assert_eq!(
            entry.tags, sketch.genome_kmers,
            "adapter keeps k-mers as-is (dedup happens downstream)"
        );
        assert_eq!(
            entry.tag_lengths,
            vec![31u8; 6],
            "every tag length must equal k"
        );
        assert_eq!(entry.sequence_id, "genome1");
        assert_eq!(entry.genome_source, "genome1");
        assert!(entry.tag_uniqueness.is_none());
        assert!(entry.tag_seqs.is_none(), "sketches store no tag sequences");
        assert_eq!(entry.enzyme, "kmer_c200_k31");

        // The same rule on the sample side makes check_enzyme_compat pass natively.
        assert!(check_enzyme_compat(
            Some(&entry.enzyme),
            Some(&sketch_enzyme_tag(200, 31)),
            false,
            "t"
        )
        .is_ok());
        // ... and hard-errors when c/k differ between DB and sample.
        assert!(check_enzyme_compat(
            Some(&entry.enzyme),
            Some(&sketch_enzyme_tag(100, 31)),
            false,
            "t"
        )
        .is_err());

        // After the standard per-genome aggregation, duplicate k-mers are deduplicated.
        let agg = aggregate_db_by_genome(vec![entry]);
        assert_eq!(agg.len(), 1);
        let mut tags = agg[0].tags.clone();
        tags.sort_unstable();
        assert_eq!(tags, vec![1, 5, 9, 42]);
        assert_eq!(agg[0].tag_lengths, vec![31u8; 4]);
        assert_eq!(agg[0].enzyme, "kmer_c200_k31");
    }

    #[test]
    fn test_sequences_sketch_view() {
        let mut kmer_counts: FxHashMap<Hash, u32> = FxHashMap::default();
        kmer_counts.insert(11, 3);
        kmer_counts.insert(22, 7);
        let sketch = SequencesSketch {
            kmer_counts,
            file_name: "sampleA.fq".to_string(),
            c: 200,
            k: 31,
            paired: false,
            sample_name: Some("sampleA".to_string()),
            mean_read_length: 150.0,
        };
        let (source, counts, total) = sketch_sample_view(&sketch);
        assert_eq!(source, "sampleA");
        assert_eq!(counts.len(), 2);
        assert_eq!(counts.get(&11), Some(&3));
        assert_eq!(
            total, 10,
            "total observations = sum of k-mer multiplicities"
        );

        // sample_name falls back to file_name when absent.
        let mut unnamed = sketch.clone();
        unnamed.sample_name = None;
        assert_eq!(sketch_sample_view(&unnamed).0, "sampleA.fq");
    }

    #[test]
    fn test_joint_ani_single_class_equals_sylph_pow() {
        // With a single length class ℓ the information-weighted least squares degenerates
        // to adj^(1/ℓ) — sylph's ANI = containment_adj^(1/k) for k-mer sketches (ℓ = k = 31).
        for &(n, d) in &[(1000usize, 1000usize), (1000, 900), (500, 250), (2000, 3)] {
            let mut classes: FxHashMap<u8, (usize, usize)> = FxHashMap::default();
            classes.insert(31u8, (n, d));
            for lambda in [None, Some(0.5), Some(2.0)] {
                let ani = joint_ani_from_classes(&classes, lambda, 0, 0.0)
                    .expect("single non-empty class must yield an ANI");
                let c = d as f64 / n as f64;
                let c_adj = match lambda {
                    Some(lam) => (c / (1.0 - (-lam).exp())).min(1.0),
                    None => c,
                };
                let expected = c_adj.powf(1.0 / 31.0);
                assert!(
                    (ani - expected).abs() < 1e-12,
                    "n={} d={} lambda={:?}: joint={} expected adj^(1/31)={}",
                    n,
                    d,
                    lambda,
                    ani,
                    expected
                );
            }
        }
    }

    #[test]
    fn test_sketch_mode_rejects_mismatch() {
        // --sketch-mode with --mismatch 1 must fail fast (before any file IO), because
        // k-mer sketches do not store tag sequences needed for neighbor generation.
        fn make_args(mismatch: usize) -> ProfileArgs {
            ProfileArgs {
                sample_file: "nonexistent.sp".to_string(),
                db_file: "nonexistent.db".to_string(),
                minimum_ani: None,
                threads: 1,
                mismatch,
                read_error_rate: None,
                no_error_correction: false,
                allow_enzyme_mismatch: false,
                out_file_name: None,
                log_path: None,
                tsv_name: "abundance_matrix.tsv".to_string(),
                taxonomy_file: None,
                gscore_threshold: 10.0,
                min_eff_coverage: 0.0,
                min_shared_tags: MIN_SHARED_TAGS_FOR_PROFILE,
                min_tags_genome: MIN_TAGS_FOR_GENOME,
                reassign_protection: false,
                sketch_mode: true,
            }
        }
        let err = profile(make_args(1)).expect_err("sketch mode + mismatch 1 must be a hard error");
        assert!(
            err.to_string().contains("--mismatch 0"),
            "error should explain the mismatch restriction, got: {}",
            err
        );

        // mismatch 0 passes the guard and fails later on the missing files instead.
        let err_ok = profile(make_args(0)).expect_err("missing files still error");
        assert!(
            !err_ok.to_string().contains("--mismatch 0"),
            "mismatch guard must not fire for --mismatch 0, got: {}",
            err_ok
        );
    }

    fn dummy_result(shared: usize, ref_tags: usize, ani: f64) -> QueryResult {
        QueryResult {
            sample_file: String::new(),
            genome_file: String::new(),
            adjusted_ani: ani,
            eff_cov: 0.0,
            ani_percentile: (0.0, 0.0),
            eff_lambda: 0.0,
            lambda_percentile: (0.0, 0.0),
            median_cov: 0.0,
            mean_cov_geq1: 0.0,
            containment_ind: String::new(),
            naive_ani: 0.0,
            contig_name: String::new(),
            ref_tags,
            shared_tags: shared,
            query_tags: 0,
            taxonomic_abundance: 0.0,
            sequence_abundance: 0.0,
        }
    }

    #[test]
    fn test_min_shared_tag_guard() {
        // A result with enough shared tags and high enough ANI should pass.
        let ok_result = dummy_result(MIN_SHARED_TAGS_FOR_PROFILE, MIN_TAGS_FOR_GENOME, PROFILE_MIN_ANI + 1.0);
        assert!(filter_results_for_profile(&ok_result, None, MIN_SHARED_TAGS_FOR_PROFILE, MIN_TAGS_FOR_GENOME));

        // A result with too few shared tags should be filtered out, even if ANI is high.
        let low_tag_result = dummy_result(MIN_SHARED_TAGS_FOR_PROFILE - 1, MIN_TAGS_FOR_GENOME, PROFILE_MIN_ANI + 1.0);
        assert!(!filter_results_for_profile(&low_tag_result, None, MIN_SHARED_TAGS_FOR_PROFILE, MIN_TAGS_FOR_GENOME));

        // Lowering the threshold should admit the same low-tag result.
        assert!(filter_results_for_profile(&low_tag_result, None, 1, MIN_TAGS_FOR_GENOME));

        // A small genome (few reference tags) is filtered at the default threshold but
        // admitted when min_tags_genome is lowered.
        let small_genome_result = dummy_result(MIN_SHARED_TAGS_FOR_PROFILE, MIN_TAGS_FOR_GENOME - 1, PROFILE_MIN_ANI + 1.0);
        assert!(!filter_results_for_profile(&small_genome_result, None, MIN_SHARED_TAGS_FOR_PROFILE, MIN_TAGS_FOR_GENOME));
        assert!(filter_results_for_profile(&small_genome_result, None, MIN_SHARED_TAGS_FOR_PROFILE, 20));

        // A result with zero shared tags should always be filtered out.
        let zero_result = dummy_result(0, MIN_TAGS_FOR_GENOME, PROFILE_MIN_ANI + 1.0);
        assert!(!filter_results_for_profile(&zero_result, None, 1, MIN_TAGS_FOR_GENOME));
    }

    #[test]
    fn test_lookup_tag_coverage_exact_and_mismatch() {
        let mut sample_counts: FxHashMap<Hash, u32> = FxHashMap::default();
        let ref_seq = b"AAACCC";
        let exact_hash = hash_bytes(ref_seq);
        let neighbor_seq = b"AAACCG"; // 1 mismatch
        let neighbor_hash = hash_bytes(neighbor_seq);

        sample_counts.insert(exact_hash, 5);
        sample_counts.insert(neighbor_hash, 3);

        // exact mode：只命中 exact
        assert_eq!(lookup_tag_coverage(exact_hash, Some(&ref_seq.to_vec()), &sample_counts, 0), Some(5));
        // mismatch mode：优先 exact（5 > 3）
        assert_eq!(lookup_tag_coverage(exact_hash, Some(&ref_seq.to_vec()), &sample_counts, 1), Some(5));

        // 把 exact 计数设为 0，mismatch mode 应命中 neighbor
        sample_counts.insert(exact_hash, 0);
        assert_eq!(lookup_tag_coverage(exact_hash, Some(&ref_seq.to_vec()), &sample_counts, 1), Some(3));

        // 无 tag_seqs 时 mismatch mode 无法生成邻居，只能 exact match
        assert_eq!(lookup_tag_coverage(exact_hash, None, &sample_counts, 1), None);
    }

    #[test]
    fn test_calculate_statistics_mismatch_increases_shared() {
        // 构造 100 个长度 32 的参考 tag，其中 70 个在样本中被观察到。
        // 在 mismatch 模型下，同样的 observed fraction 对应更低的真实 ANI。
        let total_ref = 100usize;
        let covs = vec![1u32; 70];
        let tag_lengths = vec![32u8; total_ref];

        let exact_result = calculate_statistics(covs.clone(), &tag_lengths, 1000, total_ref, 0, 0.0);
        let mm_result = calculate_statistics(covs, &tag_lengths, 1000, total_ref, 1, 0.0);

        assert_eq!(exact_result.shared_tags, mm_result.shared_tags);
        assert!(exact_result.adjusted_ani > mm_result.adjusted_ani,
            "exact ANI {} should be > mismatch ANI {}", exact_result.adjusted_ani, mm_result.adjusted_ani);
    }

    #[test]
    fn test_error_aware_inversion_recovers_true_ani() {
        // 合成场景：真实 ANI a=0.97，per-base 错误率 e=0.005，tag 长 ℓ=32。
        // 观察模型：q = a(1-e)，C = P(≤1 mismatch | q)。error-aware 反推应恢复 a。
        let a_true = 0.97;
        let e = 0.005;
        let len = 32usize;
        let q = a_true * (1.0 - e);
        let c = crate::extract::p_detect_one_mismatch(q, len);

        let a_recovered = ani_from_containment_one_mismatch_err(c, len, e);
        assert!((a_recovered - a_true).abs() < 0.002,
            "error-aware inversion recovered {}, expected ~{}", a_recovered, a_true);

        // 未校正（旧）反推应明显偏低（约 (1-e) 因子）。
        let a_uncorrected = ani_from_containment_one_mismatch(c, len);
        assert!(a_uncorrected < a_true - 0.003,
            "uncorrected inversion {} should be biased well below {}", a_uncorrected, a_true);
    }

    #[test]
    fn test_error_aware_inversion_e0_matches_old() {
        // e=0 时必须与旧公式逐位一致。
        for &len in &[30usize, 32, 34] {
            for &c in &[0.05, 0.3, 0.7, 0.95, 1.0] {
                let old = ani_from_containment_one_mismatch(c, len);
                let new = ani_from_containment_one_mismatch_err(c, len, 0.0);
                assert!((old - new).abs() < 1e-12, "len={} c={}: {} vs {}", len, c, old, new);
            }
        }
    }

    #[test]
    fn test_pooled_error_rate_recovers_e() {
        // ℓ=32，q=0.995（e=0.005）：E[mm1-only]/E[exact] = ℓ(1-q)/q。
        let len = 32u8;
        let q = 0.995f64;
        let exact = 10_000u64;
        let mm1 = (len as f64 * (1.0 - q) / q * exact as f64).round() as u64;
        let mut classes: FxHashMap<u8, (u64, u64)> = FxHashMap::default();
        classes.insert(len, (exact, mm1));
        let e = pooled_error_rate(&classes).expect("should estimate");
        assert!((e - 0.005).abs() < 0.0005, "estimated e={}", e);

        // tag 太少 → None（调用方回退 e=0）。
        let mut sparse: FxHashMap<u8, (u64, u64)> = FxHashMap::default();
        sparse.insert(len, (10, 1));
        assert!(pooled_error_rate(&sparse).is_none());
    }

    #[test]
    fn test_sample_error_rate_override_paths() {
        // --no-error-correction 恒为 0（即使有 override）；--read-error-rate 优先；
        // 无 override 且 mm1 时走数据估计。这里用空 DB/样本触发估计回退（e=0）。
        let empty_db: Vec<SyldbEntry> = Vec::new();
        let empty_counts: FxHashMap<Hash, u32> = FxHashMap::default();

        assert_eq!(sample_error_rate(&empty_db, &empty_counts, 1, Some(0.01), true), 0.0);
        assert_eq!(sample_error_rate(&empty_db, &empty_counts, 1, Some(0.01), false), 0.01);
        // mismatch=0 时完全忽略 error correction 设置。
        assert_eq!(sample_error_rate(&empty_db, &empty_counts, 0, Some(0.01), false), 0.0);
        // 无 override + 无数据 → 估计回退到 0（旧行为）。
        assert_eq!(sample_error_rate(&empty_db, &empty_counts, 1, None, false), 0.0);
    }

    #[test]
    fn test_estimate_read_error_rate_synthetic() {
        // 构造一个"近同一"基因组：1000 个长度 32 的参考 tag。
        // 样本中 exact 命中 900 个，另有 144 个以 1-mismatch 邻居形式命中
        // （比例 ≈ ℓ(1-q)/q = 32·0.005/0.995 ≈ 0.1608 → e 应估计为 ~0.005）。
        let n_tags = 1000usize;
        let n_exact = 900usize;
        let n_mm1 = 144usize;
        let mut tags = Vec::with_capacity(n_tags);
        let mut tag_seqs = Vec::with_capacity(n_tags);
        let mut sample_counts: FxHashMap<Hash, u32> = FxHashMap::default();

        // 与 extract 侧一致的 canonical 形式（正向与反向互补取字典序较小者）。
        let revcomp = |s: &[u8]| -> Vec<u8> {
            s.iter().rev().map(|b| match b { b'A' => b'T', b'T' => b'A', b'C' => b'G', _ => b'C' }).collect()
        };
        let canon = |s: &[u8]| -> Vec<u8> {
            let rc = revcomp(s);
            if s <= rc.as_slice() { s.to_vec() } else { rc }
        };

        // 用确定性伪随机序列作为 tag（互不相同的合法 ACGT 序列即可）。
        let mut rng: u64 = 0x9e3779b97f4a7c15;
        let mut mk_seq = move || -> Vec<u8> {
            (0..32).map(|_| {
                rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGT"[(rng >> 33) as usize & 3]
            }).collect()
        };

        for _ in 0..n_exact {
            let seq = canon(&mk_seq());
            tags.push(hash_bytes(&seq));
            sample_counts.insert(hash_bytes(&seq), 1);
            tag_seqs.push(seq);
        }
        // mm1-only：把序列的第 5 位换成不同碱基，样本里只存突变体的 canonical hash。
        for _ in 0..n_mm1 {
            let seq = canon(&mk_seq());
            let mut mutated = seq.clone();
            mutated[5] = match mutated[5] { b'A' => b'C', b'C' => b'G', b'G' => b'T', _ => b'A' };
            let m_hash = hash_bytes(&canon(&mutated));
            // one_mismatch_canonical_hashes 应覆盖该突变体；防御性检查，不满足就跳过。
            if one_mismatch_canonical_hashes(&seq).contains(&m_hash) {
                tags.push(hash_bytes(&seq));
                sample_counts.insert(m_hash, 1);
                tag_seqs.push(seq);
            }
        }
        // 其余 tag 不在样本中。
        while tags.len() < n_tags {
            let seq = canon(&mk_seq());
            tags.push(hash_bytes(&seq));
            tag_seqs.push(seq);
        }

        let entry = SyldbEntry {
            sequence_id: "g1".to_string(),
            tags,
            tag_lengths: vec![32u8; n_tags],
            genome_source: "g1".to_string(),
            tag_uniqueness: None,
            tag_seqs: Some(tag_seqs),
            enzyme: "BcgI".to_string(),
        };
        let db = vec![entry];
        let e = estimate_read_error_rate(&db, &sample_counts);
        assert!((e - 0.005).abs() < 0.002, "estimated e={}", e);
    }

    #[test]
    fn test_enzyme_compat_guard() {
        // 匹配（含顺序/空白差异）→ 通过。
        assert!(check_enzyme_compat(Some("BcgI"), Some("BcgI"), false, "t").is_ok());
        assert!(check_enzyme_compat(Some("BcgI,BslFI"), Some(" BslFI , BcgI "), false, "t").is_ok());

        // DB 有样本缺失的酶 → 硬错误；--allow-enzyme-mismatch → 通过。
        assert!(check_enzyme_compat(Some("BcgI,BslFI,CspCI,AloI"), Some("BcgI"), false, "t").is_err());
        assert!(check_enzyme_compat(Some("BcgI,BslFI,CspCI,AloI"), Some("BcgI"), true, "t").is_ok());

        // 样本多出 DB 没有的酶 → 仅警告，不硬错。
        assert!(check_enzyme_compat(Some("BcgI"), Some("BcgI,AloI"), false, "t").is_ok());

        // 任一侧缺酶信息（旧格式文件）→ 跳过检查。
        assert!(check_enzyme_compat(None, Some("BcgI"), false, "t").is_ok());
        assert!(check_enzyme_compat(Some("BcgI"), None, false, "t").is_ok());
        assert!(check_enzyme_compat(Some(""), Some("BcgI"), false, "t").is_ok());
    }

    /// 构造过度剥离场景（v2 阈值）：
    /// - A = 赢家，与 B 共享 tag 1..=100、与 C/D 共享 1..=60，ANI 最高 → 赢得全部共享 tag。
    /// - B = 近全灭输家：100 个共享 tag 全丢，仅剩 4 个独有 tag（存活 4% ≤ 5%）→ 应被保护。
    /// - C = 普通菌株簇输家：丢 60/100（存活 40% > 5%）→ 不应被保护（v1 的 >50% 丢失
    ///   阈值会误保护它）。
    /// - D = 弱初始证据：60 个共享 tag 全丢（存活 0%），但初始共享 60 < 100 → 不应被保护。
    fn over_stripped_fixture() -> (Vec<SyldbEntry>, Vec<SylspEntry>, Vec<QueryResult>) {
        let shared_ab: Vec<Hash> = (1u64..=100).collect();
        let b_unique: Vec<Hash> = (2001u64..=2004).collect();
        let c_unique: Vec<Hash> = (3001u64..=3040).collect();
        let a_unique: Hash = 9001; // 不在样本中
        let d_unique: Hash = 9501; // 不在样本中

        let mk_db = |id: &str, shared: &[Hash], unique: &[Hash]| {
            let mut tags: Vec<Hash> = shared.to_vec();
            tags.extend_from_slice(unique);
            SyldbEntry {
                sequence_id: id.to_string(),
                tag_lengths: vec![32u8; tags.len()],
                tags,
                genome_source: format!("{}.fasta", id),
                tag_uniqueness: None,
                tag_seqs: None,
                enzyme: "BcgI".to_string(),
            }
        };
        let db = vec![
            mk_db("A", &shared_ab, &[a_unique]),
            mk_db("B", &shared_ab, &b_unique),
            mk_db("C", &shared_ab[..60], &c_unique),
            mk_db("D", &shared_ab[..60], &[d_unique]),
        ];

        // 样本：共享 tag 覆盖度 4，B 独有 tag 覆盖度 5，C 独有 tag 覆盖度 3。
        let mut sample = Vec::new();
        let mut push_tags = |tags: &[Hash], copies: usize, prefix: &str| {
            for &t in tags {
                for i in 0..copies {
                    sample.push(SylspEntry {
                        sequence_id: format!("{}{}_{}", prefix, t, i),
                        tag: t,
                        tag_length: 32,
                        sample_source: "s1".to_string(),
                        enzyme: "BcgI".to_string(),
                    });
                }
            }
        };
        push_tags(&shared_ab, 4, "rs");
        push_tags(&b_unique, 5, "rb");
        push_tags(&c_unique, 3, "rc");

        // 初始（重分配前）结果：均通过 profile 过滤，A 的 ANI 最高 → 赢得所有共享 tag。
        let mk_init = |id: &str, shared: usize, ref_tags: usize, ani: f64| {
            let mut r = dummy_result(shared, ref_tags, ani);
            r.contig_name = id.to_string();
            r.sample_file = "s1".to_string();
            r
        };
        let initial = vec![
            mk_init("A", 100, 101, 99.9),
            mk_init("B", 100, 104, 99.0),
            mk_init("C", 100, 100, 98.5),
            mk_init("D", 60, 61, 98.0),
        ];

        (db, sample, initial)
    }

    #[test]
    fn test_reassign_over_stripping_protection() {
        let (db, sample, initial) = over_stripped_fixture();
        let (winner_map, interner) = build_winner_table(&initial, &db, false);

        let min_ani = 95.0;
        let min_shared = 50;
        let min_tags_genome = 1;

        // 保护开启：B（存活 4%）被保留 —— 初始 ANI/shared_tags，重分配后的 eff_cov。
        let protected = recalculate_with_winner_table(
            &db,
            &sample,
            &winner_map,
            &interner,
            min_ani,
            false,
            0,
            min_shared,
            min_tags_genome,
            None,
            false,
            &initial,
            true,
        );
        let b = protected
            .iter()
            .find(|r| r.genome_file == "B.fasta")
            .expect("near-annihilated genome B (4% survival) should be protected");
        assert_eq!(
            b.adjusted_ani, 99.0,
            "protected genome keeps pre-reassignment ANI"
        );
        assert_eq!(
            b.shared_tags, 100,
            "protected genome keeps pre-reassignment shared tags"
        );
        assert_eq!(
            b.eff_cov, 5.0,
            "protected genome uses post-reassignment coverage for abundance"
        );
        assert!(
            protected.iter().any(|r| r.genome_file == "A.fasta"),
            "winner A still reported"
        );
        assert!(
            !protected.iter().any(|r| r.genome_file == "C.fasta"),
            "ordinary cluster loser C (60% loss, 40% survival) must NOT be protected"
        );
        assert!(
            !protected.iter().any(|r| r.genome_file == "D.fasta"),
            "genome D with weak initial evidence (60 < 100 initial shared) must NOT be protected"
        );

        // 保护关闭（--no-reassign-protection）：B/C/D 均被剥离后过滤掉，与旧行为一致。
        let unprotected = recalculate_with_winner_table(
            &db,
            &sample,
            &winner_map,
            &interner,
            min_ani,
            false,
            0,
            min_shared,
            min_tags_genome,
            None,
            false,
            &initial,
            false,
        );
        assert!(
            !unprotected.iter().any(|r| r.genome_file == "B.fasta"),
            "without protection B is stripped below the floor and dropped"
        );
        assert!(!unprotected.iter().any(|r| r.genome_file == "C.fasta"));
        assert!(!unprotected.iter().any(|r| r.genome_file == "D.fasta"));
        assert!(unprotected.iter().any(|r| r.genome_file == "A.fasta"));

        // 防膨胀护栏：初始 ANI 低于 min_ani 的基因组即使被近全灭也不受保护。
        let mut low_ani_initial = initial.clone();
        // B 初始 ANI 不达标；重建 winner table（B 的 ANI 变了，但 A 仍赢所有共享 tag）。
        low_ani_initial[1].adjusted_ani = 90.0;
        let (winner_map2, interner2) = build_winner_table(&low_ani_initial, &db, false);
        let guarded = recalculate_with_winner_table(
            &db,
            &sample,
            &winner_map2,
            &interner2,
            min_ani,
            false,
            0,
            min_shared,
            min_tags_genome,
            None,
            false,
            &low_ani_initial,
            true,
        );
        assert!(
            !guarded.iter().any(|r| r.genome_file == "B.fasta"),
            "genome failing min-ani on pre-reassignment stats must not be protected"
        );
    }
}
