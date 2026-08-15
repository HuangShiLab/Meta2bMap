use anyhow::{Context, Result};
use bio::io::{fasta, fastq};
use needletail::parse_fastx_file;
use flate2::{read::GzDecoder, write::GzEncoder, Compression};
// bytes::Regex 直接在 &[u8] 上匹配，避免每条序列的 from_utf8_lossy 全量拷贝。
// 所有酶切模式都是 ASCII 的 [ACGT] 字符类，字节匹配与字符匹配完全等价。
use regex::bytes::Regex;
use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Read, Write},
    path::{Path, PathBuf},
    thread,
    time::Duration,
};

use crate::cmdline::ExtractArgs;
use serde::{Serialize, Deserialize};
use rayon::prelude::*;
use std::sync::{Arc, Mutex};
use crate::constants::{Hash, hash_bytes};
// 添加fxhash导入
use fxhash::{FxHashMap, FxHashSet};

// 添加内存统计导入
use memory_stats::memory_stats;
use log::*;

// 类型别名，与sylph保持一致
pub type TagHash = Vec<u8>;
pub type SampleId = String;
pub type SampleStatsMap = FxHashMap<SampleId, ExtractionStats>;

// 优化的压缩设置
fn get_optimal_compression() -> Compression {
    // 根据系统性能调整压缩级别
    if std::env::var("FAST_COMPRESSION").is_ok() {
        Compression::fast()
    } else if std::env::var("BEST_COMPRESSION").is_ok() {
        Compression::best()
    } else {
        Compression::default()
    }
}

// 优化的文件大小检测
fn get_file_size_optimized(path: &Path) -> Result<u64> {
    let metadata = std::fs::metadata(path)?;
    Ok(metadata.len())
}

// 优化的缓冲区大小计算
fn calculate_optimal_buffer_size(file_size: u64, is_compressed: bool) -> usize {
    let base_size = if file_size > 1024 * 1024 * 1024 {
        // 大文件：256KB
        256 * 1024
    } else if file_size > 100 * 1024 * 1024 {
        // 中等文件：128KB
        128 * 1024
    } else {
        // 小文件：64KB
        64 * 1024
    };
    
    if is_compressed {
        base_size * 2
    } else {
        base_size
    }
}



// 内存监控函数，参考sketch中的check_vram_and_block
// 使用 physical_mem 而非 virtual_mem，避免 macOS 上虚拟地址空间过大导致死锁。
// 注意：阻塞期间所有并行线程可能互相等待（谁都不释放内存）形成死锁，
// 因此阻塞必须有上限和可见的告警——宁可超时后继续（最坏是 OOM，错误可见），
// 也不能无声地无限挂起。
pub fn check_vram_and_block(max_ram: usize, file: &str) {
    if let Some(usage) = memory_stats() {
        let mut gb_usage_curr = usage.physical_mem as f64 / 1_000_000_000 as f64;
        if (max_ram as f64) < gb_usage_curr {
            eprintln!(
                "WARNING: memory limit reached ({:.1} GB > {} GB). Blocking extract for {} until memory frees. If this persists, re-run with a higher --max-ram.",
                gb_usage_curr, max_ram, file
            );
        }
        let mut blocked_secs = 0u64;
        while (max_ram as f64) < gb_usage_curr {
            thread::sleep(Duration::from_secs(5));
            blocked_secs += 5;
            if let Some(usage) = memory_stats() {
                gb_usage_curr = usage.physical_mem as f64 / 1_000_000_000 as f64;
            } else {
                break;
            }
            if blocked_secs >= 600 {
                eprintln!(
                    "WARNING: extract for {} blocked for >10 min at {:.1} GB (limit {} GB); worker threads may be waiting on each other (deadlock). Proceeding anyway; consider a higher --max-ram.",
                    file, gb_usage_curr, max_ram
                );
                break;
            }
            if blocked_secs % 60 == 0 {
                eprintln!(
                    "WARNING: extract for {} still blocked: {:.1} GB > {} GB limit",
                    file, gb_usage_curr, max_ram
                );
            }
        }
    }
}

// 默认内存上限：总物理内存的 75%。旧版硬编码 16GB 在现代机器上会频繁误触发
// 上面的阻塞 guard——多个并行线程同时阻塞会互相等待形成死锁（HPC 建库曾因此
// 无声挂死数小时）。下限 7GB 与 extract() 的显式校验保持一致。
pub fn default_max_ram_gb() -> usize {
    use sysinfo::SystemExt;
    let mut sys = sysinfo::System::new();
    sys.refresh_memory();
    let total_gb = sys.total_memory() as f64 / 1_000_000_000.0;
    ((total_gb * 0.75) as usize).max(7)
}

// 动态内存管理函数
pub fn get_memory_usage() -> Option<f64> {
    memory_stats().map(|usage| usage.physical_mem as f64 / 1_000_000_000 as f64)
}







// 内存安全的处理函数
pub fn safe_process_with_memory_check<F, T>(
    max_ram: usize,
    file: &str,
    process_fn: F,
) -> Result<T>
where
    F: FnOnce() -> Result<T>,
{
    // 检查当前内存使用
    if let Some(current_memory) = get_memory_usage() {
        if current_memory > max_ram as f64 {
            check_vram_and_block(max_ram, file);
        }
    }
    
    // 执行处理函数
    let result = process_fn()?;
    
    // 处理完成后再次检查内存
    if let Some(current_memory) = get_memory_usage() {
        if current_memory > max_ram as f64 * 0.8 {
            // 如果内存使用超过80%，等待一下让系统回收内存
            thread::sleep(Duration::from_millis(100));
        }
    }
    
    Ok(result)
}

pub const ENZYME_DEFINITIONS: &[(&str, &[&str])] = &[
    ("CspCI", &[
        r"[ACGT]{11}CAA[ACGT]{5}GTGG[ACGT]{10}",
        r"[ACGT]{10}CCAC[ACGT]{5}TTG[ACGT]{11}",
    ]),
    ("AloI", &[
        r"[ACGT]{7}GAAC[ACGT]{6}TCC[ACGT]{7}",
        r"[ACGT]{7}GGA[ACGT]{6}GTTC[ACGT]{7}",
    ]),
    ("BsaXI", &[
        r"[ACGT]{9}AC[ACGT]{5}CTCC[ACGT]{7}",
        r"[ACGT]{7}GGAG[ACGT]{5}GT[ACGT]{9}",
    ]),
    ("BaeI", &[
        r"[ACGT]{10}AC[ACGT]{4}GTA[CT]C[ACGT]{7}",
        r"[ACGT]{7}G[AG]TAC[ACGT]{4}GT[ACGT]{10}",
    ]),
    ("BcgI", &[
        r"[ACGT]{10}CGA[ACGT]{6}TGC[ACGT]{10}",
        r"[ACGT]{10}GCA[ACGT]{6}TCG[ACGT]{10}",
    ]),
    ("CjeI", &[
        r"[ACGT]{8}CCA[ACGT]{6}GT[ACGT]{9}",
        r"[ACGT]{9}AC[ACGT]{6}TGG[ACGT]{8}",
    ]),
    ("PpiI", &[
        r"[ACGT]{7}GAAC[ACGT]{5}CTC[ACGT]{8}",
        r"[ACGT]{8}GAG[ACGT]{5}GTTC[ACGT]{7}",
    ]),
    ("PsrI", &[
        r"[ACGT]{7}GAAC[ACGT]{6}TAC[ACGT]{7}",
        r"[ACGT]{7}GTA[ACGT]{6}GTTC[ACGT]{7}",
    ]),
    ("BplI", &[
        r"[ACGT]{8}GAG[ACGT]{5}CTC[ACGT]{8}",
    ]),
    ("FalI", &[
        r"[ACGT]{8}AAG[ACGT]{5}CTT[ACGT]{8}",
    ]),
    ("Bsp24I", &[
        r"[ACGT]{8}GAC[ACGT]{6}TGG[ACGT]{7}",
        r"[ACGT]{7}CCA[ACGT]{6}GTC[ACGT]{8}",
    ]),
    ("HaeIV", &[
        r"[ACGT]{7}GA[CT][ACGT]{5}[AG]TC[ACGT]{9}",
        r"[ACGT]{9}GA[CT][ACGT]{5}[AG]TC[ACGT]{7}",
    ]),
    ("CjePI", &[
        r"[ACGT]{7}CCA[ACGT]{7}TC[ACGT]{8}",
        r"[ACGT]{8}GA[ACGT]{7}TGG[ACGT]{7}",
    ]),
    ("Hin4I", &[
        r"[ACGT]{8}GA[CT][ACGT]{5}[GAC]TC[ACGT]{8}",
        r"[ACGT]{8}GA[CTG][ACGT]{5}[AG]TC[ACGT]{8}",
    ]),
    ("AlfI", &[
        r"[ACGT]{10}GCA[ACGT]{6}TGC[ACGT]{10}",
    ]),
    ("BslFI", &[
        r"[ACGT]{6}GGGAC[ACGT]{14}",
        r"[ACGT]{14}GTCCC[ACGT]{6}",
    ]),
];

// 定义每个内切酶的标签长度（固定匹配碱基数 + 自由匹配碱基数）
// 数值等于该酶 @site 正则实际匹配的总长度（与 2bRADExtraction.pl 完全一致）。
// 注意：AloI / BaeI / HaeIV / Hin4I 原先的数值与其自身正则模式的真实匹配长度不符
// （旧值分别为 20/27/25/25，均属手误漏算），导致下方 extract_and_validate_tags /
// extract_tags_avx2 里 "matched.len() > tag_length 时取中间窗口" 的逻辑会把命中的
// 完整酶切位点错误地截短，与 Perl 版本输出的 tag 不一致，现已修正。
pub const ENZYME_TAG_LENGTHS: &[(&str, usize)] = &[
    ("CspCI", 33),  // 11 + 3 + 5 + 4 + 10 = 33
    ("AloI", 27),   // 7 + 4 + 6 + 3 + 7 = 27
    ("BsaXI", 27),  // 9 + 2 + 5 + 4 + 7 = 27
    ("BaeI", 28),   // 10 + 2 + 4 + 3 + 1 + 1 + 7 = 28 (GTA[CT]C = 5 literal/degenerate positions)
    ("BcgI", 32),   // 10 + 3 + 6 + 3 + 10 = 32
    ("CjeI", 28),   // 8 + 3 + 6 + 2 + 9 = 28
    ("PpiI", 27),   // 7 + 4 + 5 + 3 + 8 = 27
    ("PsrI", 27),   // 7 + 4 + 6 + 3 + 7 = 27
    ("BplI", 27),   // 8 + 3 + 5 + 3 + 8 = 27
    ("FalI", 27),   // 8 + 3 + 5 + 3 + 8 = 27
    ("Bsp24I", 27), // 8 + 3 + 6 + 3 + 7 = 27
    ("HaeIV", 27),  // 7 + 3 + 5 + 3 + 9 = 27 (GA[CT] = 3, [AG]TC = 3)
    ("CjePI", 27),  // 7 + 3 + 7 + 2 + 8 = 27
    ("Hin4I", 27),  // 8 + 3 + 5 + 3 + 8 = 27 (GA[CT] = 3, [GAC]TC = 3)
    ("AlfI", 32),   // 10 + 3 + 6 + 3 + 10 = 32
    ("BslFI", 25),  // 6 + 5 + 14 = 25
];

#[derive(Debug)]
pub struct EnzymeSpec {
    pub name: String,
    pub patterns: Vec<Regex>,
    /// 每个 pattern 对应的 tag 长度。多酶模式下 pattern 来自不同酶，长度可能不同。
    pub pattern_tag_lengths: Vec<usize>,
    /// 主 tag 长度（日志/统计使用，取第一个酶的 tag 长度）。
    pub tag_length: usize,
}

impl EnzymeSpec {
    pub fn new(name: &str) -> Result<Self> {
        let names: Vec<&str> = if name.eq_ignore_ascii_case("all") {
            ENZYME_DEFINITIONS.iter().map(|(n, _)| *n).collect()
        } else {
            name.split(',').map(|s| s.trim()).filter(|s| !s.is_empty()).collect()
        };

        if names.is_empty() {
            return Err(anyhow::anyhow!("No enzyme specified"));
        }

        let mut all_patterns = Vec::new();
        let mut all_lengths = Vec::new();
        let mut seen = FxHashSet::default();

        for enzyme_name in &names {
            if !seen.insert(*enzyme_name) {
                continue;
            }
            let def = ENZYME_DEFINITIONS
                .iter()
                .find(|(e, _)| *e == *enzyme_name)
                .ok_or_else(|| anyhow::anyhow!("Unsupported enzyme: {}", enzyme_name))?;

            let tag_length = ENZYME_TAG_LENGTHS
                .iter()
                .find(|(n, _)| *n == def.0)
                .map(|(_, len)| *len)
                .ok_or_else(|| anyhow::anyhow!("Missing tag length for enzyme: {}", def.0))?;

            for pat in def.1 {
                all_patterns.push(Regex::new(pat).context(format!("Invalid regex pattern: {}", pat))?);
                all_lengths.push(tag_length);
            }
        }

        let primary_length = all_lengths[0];
        let display_name = names.join(",");

        Ok(Self {
            name: display_name,
            patterns: all_patterns,
            pattern_tag_lengths: all_lengths,
            tag_length: primary_length,
        })
    }
}

/// tag 最大长度（当前最长的 CspCI = 33），canonical 化时用作栈上缓冲区大小。
const MAX_TAG_LEN: usize = 40;

/// 碱基互补查表，替代逐碱基 match。非 ACGT 原样保留。
const COMPLEMENT: [u8; 256] = {
    let mut table = [0u8; 256];
    let mut i = 0;
    while i < 256 {
        table[i] = i as u8;
        i += 1;
    }
    table[b'A' as usize] = b'T';
    table[b'T' as usize] = b'A';
    table[b'C' as usize] = b'G';
    table[b'G' as usize] = b'C';
    table
};

/// 将 tag 的 canonical 形式（正向与反向互补中字典序较小者）写入 `buf`，返回长度。
/// 全程零堆分配：反向互补写入 buf，若正向更小则用正向覆盖。
#[inline]
fn canonicalize_into(tag: &[u8], buf: &mut [u8; MAX_TAG_LEN]) -> usize {
    let n = tag.len();
    for i in 0..n {
        buf[i] = COMPLEMENT[tag[n - 1 - i] as usize];
    }
    if tag <= &buf[..n] {
        buf[..n].copy_from_slice(tag);
    }
    n
}

/// 从序列中提取所有 canonical tag 的哈希及其长度（按 u64 去重），不为每个 tag 分配 Vec。
/// 下游（syldb/sylsp）只使用 tag 的 u64 哈希，因此这是最热路径的首选。
fn extract_tag_hashes(seq: &[u8], enzyme: &EnzymeSpec) -> Vec<(Hash, u8)> {
    let mut hashes = Vec::with_capacity(64);
    let mut seen = FxHashSet::default();
    let mut buf = [0u8; MAX_TAG_LEN];
    for (offset, len) in find_all_tag_positions(seq, enzyme) {
        let n = canonicalize_into(&seq[offset..offset + len], &mut buf);
        let h = hash_bytes(&buf[..n]);
        if seen.insert(h) {
            hashes.push((h, n as u8));
        }
    }
    hashes
}

/// 与 `extract_tag_hashes` 相同的扫描，但返回 canonical tag 的字节序列及其长度（按哈希去重）。
/// 仅用于需要输出实际 tag 序列（如 reads→FASTA）的少数路径。
fn extract_canonical_tags(seq: &[u8], enzyme: &EnzymeSpec) -> Vec<(TagHash, u8)> {
    let mut tags = Vec::with_capacity(64);
    let mut seen = FxHashSet::default();
    let mut buf = [0u8; MAX_TAG_LEN];
    for (offset, len) in find_all_tag_positions(seq, enzyme) {
        let n = canonicalize_into(&seq[offset..offset + len], &mut buf);
        let h = hash_bytes(&buf[..n]);
        if seen.insert(h) {
            tags.push((buf[..n].to_vec(), n as u8));
        }
    }
    tags
}

/// 生成 `tag` 的所有 canonical 1-mismatch 变体的哈希。
/// 每个位置尝试 3 个替代碱基，并对每个变体做 canonical 化（取 forward/revcomp 字典序较小者），
/// 因此与样本提取时的 canonical 化规则一致。
pub fn one_mismatch_canonical_hashes(tag: &[u8]) -> Vec<Hash> {
    let mut out = Vec::with_capacity(tag.len() * 3);
    let mut buf = [0u8; MAX_TAG_LEN];
    let mut neighbor_buf = [0u8; MAX_TAG_LEN];
    for (i, &orig) in tag.iter().enumerate() {
        for alt in [b'A', b'C', b'G', b'T'] {
            if alt == orig {
                continue;
            }
            neighbor_buf[..tag.len()].copy_from_slice(tag);
            neighbor_buf[i] = alt;
            let n = canonicalize_into(&neighbor_buf[..tag.len()], &mut buf);
            out.push(hash_bytes(&buf[..n]));
        }
    }
    out
}

/// 在 exact matching 下，检出概率为 P_cons(a, ℓ) = a^ℓ。
pub fn p_detect_exact(a: f64, len: usize) -> f64 {
    a.powi(len as i32)
}

/// 允许 ≤1 mismatch 时，检出概率为 Σ_{j=0}^{1} C(ℓ,j) (1-a)^j a^{ℓ-j}
///                         = a^ℓ + ℓ·a^{ℓ-1}·(1-a)。
pub fn p_detect_one_mismatch(a: f64, len: usize) -> f64 {
    if len == 0 {
        return 1.0;
    }
    let l = len as f64;
    let exact = a.powi(len as i32);
    let one_off = l * a.powi((len - 1) as i32) * (1.0 - a);
    (exact + one_off).min(1.0)
}

/// 从观察到的 containment（已做覆盖度校正后）反推 ANI，假设允许 ≤1 mismatch。
/// 用二分法求解 p_detect_one_mismatch(a, len) = containment。
pub fn ani_from_containment_one_mismatch(containment: f64, len: usize) -> f64 {
    let c = containment.clamp(0.0, 1.0);
    if c <= 0.0 {
        return 0.0;
    }
    if c >= 1.0 {
        return 1.0;
    }
    let mut lo = 0.0;
    let mut hi = 1.0;
    for _ in 0..60 {
        let mid = (lo + hi) / 2.0;
        let p = p_detect_one_mismatch(mid, len);
        if p < c {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    ((lo + hi) / 2.0).clamp(0.0, 1.0)
}

/// Error-aware variant of `ani_from_containment_one_mismatch`.
///
/// Model: sample reads carry a per-base sequencing error rate `e`, so a homologous
/// tag from a genome with true per-base identity `a` is observed with per-base match
/// probability `q = a·(1-e)`, not `a`. The observed (coverage-corrected) containment
/// therefore satisfies `C = P(≤1 mismatch | q) = q^ℓ + ℓ·q^{ℓ-1}·(1-q)`.
/// Inversion: solve for `q` with the existing bisection (identical formula), then
/// recover `a = q / (1-e)`, clamped to ≤ 1. With `e = 0` this reduces exactly to
/// `ani_from_containment_one_mismatch`.
pub fn ani_from_containment_one_mismatch_err(containment: f64, len: usize, e: f64) -> f64 {
    let q = ani_from_containment_one_mismatch(containment, len);
    let denom = 1.0 - e;
    if denom <= 0.0 {
        return q; // 非法的 e：不做校正，避免除零/放大
    }
    (q / denom).clamp(0.0, 1.0)
}

/// Scan `seq_str` for every occurrence — including *overlapping* ones — of
/// any of `enzyme`'s recognition patterns, and return `(start, length)` for
/// each hit.
///
/// `2bRADExtraction.pl`'s `Electronic_enzyme`/`fastq` subroutines find every
/// occurrence of a site by rewinding the regex cursor to `match_start + 1`
/// after each hit, rather than continuing from the end of the match (as
/// `Regex::find_iter` does by default). That means two recognition sites
/// that overlap by a few bases are *both* reported in Perl, but
/// `find_iter`-based scanning would silently miss the second one.
///
/// This reproduces that exact rewind technique using `Regex::find_at`,
/// which — unlike testing an anchored copy of the pattern against every
/// fixed-length window — lets the regex engine keep using its normal
/// literal/Aho-Corasick prefiltering to jump straight to the next candidate
/// position instead of re-verifying the whole pattern at every single
/// offset. In benchmarks this is ~15-17x faster than a naive per-offset
/// windowed scan, and within noise of the original (overlap-missing)
/// `find_iter` loop it replaces. All enzyme patterns use fixed-count `{n}`
/// repetitions only (no variable-length quantifiers), so every match is
/// guaranteed to have length `enzyme.tag_length` regardless of where it
/// starts — no separate anchoring or length check is needed.
fn find_all_tag_positions(seq: &[u8], enzyme: &EnzymeSpec) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    for (pattern, &tag_len) in enzyme.patterns.iter().zip(&enzyme.pattern_tag_lengths) {
        let mut start = 0usize;
        while start <= seq.len() {
            match pattern.find_at(seq, start) {
                Some(m) => {
                    let mstart = m.start();
                    out.push((mstart, tag_len));
                    start = mstart + 1; // rewind: mirrors Perl's `pos($seq) = match_start + 1`
                }
                None => break,
            }
        }
    }
    out
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct SyldbEntry {
    pub sequence_id: String,
    pub tags: Vec<Hash>,
    /// 每个 tag 的长度（bp），与 `tags` 一一对应。多酶联合 sketch 需要按长度分区估计 ANI。
    pub tag_lengths: Vec<u8>,
    pub genome_source: String,
    // 新增字段：标记每个tag是否为unique（taxa-specific）
    pub tag_uniqueness: Option<Vec<bool>>,
    /// 每个 tag 的 canonical 序列字节。用于 error-tolerant matching（≤1 mismatch）。
    /// 旧版数据库不含此字段，反序列化时为 None，此时只能做 exact matching。
    #[serde(default)]
    pub tag_seqs: Option<Vec<TagHash>>,
    /// 构建该数据库时使用的酶（或酶组合），使 .syldb 文件自描述。
    /// 旧版数据库不含此字段，反序列化时为空字符串。
    #[serde(default)]
    pub enzyme: String,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct SylspEntry {
    pub sequence_id: String,
    pub tag: Hash,
    /// 该 tag 的长度（bp），用于多酶联合 ANI 的长度分区。
    pub tag_length: u8,
    pub sample_source: String,
    /// 提取该样本 sketch 时使用的酶（或酶组合），使 .sylsp 文件自描述，
    /// 用于 profile/query 时与 DB 的酶集合做一致性检查。
    /// 注意：bincode 不是自描述格式，旧版 .sylsp（无此字段）反序列化会失败，需要重新 extract。
    #[serde(default)]
    pub enzyme: String,
}

#[derive(Serialize, Deserialize, Debug, PartialEq, Hash, PartialOrd, Eq, Ord, Default, Clone)]
pub struct GenomeSketch {
    pub file_name: String,
    pub first_contig_name: String,
    pub gn_size: usize,
    pub c: usize,
    pub k: usize,
    pub min_spacing: usize,
    pub genome_kmers: Vec<Hash>,
}

#[derive(Serialize, Deserialize, Debug, PartialEq, Hash, PartialOrd, Eq, Ord, Default, Clone)]
pub struct GenomeSketchInspect {
    pub file_name: String,
    pub genome_kmers_num: usize,
    pub first_contig_name: String,
    pub genome_size: usize,
}

impl From<GenomeSketch> for GenomeSketchInspect {
    fn from(sk: GenomeSketch) -> Self {
        GenomeSketchInspect {
            genome_kmers_num: sk.genome_kmers.len(),
            file_name: sk.file_name,
            first_contig_name: sk.first_contig_name,
            genome_size: sk.gn_size,
        }
    }
}

#[derive(Serialize, Deserialize, Debug, PartialEq, Hash, PartialOrd, Eq, Ord, Default, Clone)]
pub struct DatabaseSketch {
    pub database_file: String,
    pub c: usize,
    pub k: usize,
    pub min_spacing_parameter: usize,
    pub genome_files: Vec<GenomeSketchInspect>,
}

#[allow(dead_code)]
pub fn process_input(
    input_files: Vec<PathBuf>,
    sample_output_dir: &Path,
    enzyme_name: &str,
    _threads: usize,
    format: &str,
) -> Result<()> {
    let enzyme = EnzymeSpec::new(enzyme_name)
        .context(format!("Unsupported enzyme: {}", enzyme_name))?;

    for input_path in &input_files {
        // 确定输入文件类型
        let is_fasta = is_fasta_file(input_path)
            .context("Failed to determine if file is FASTA")?;
        let is_fastq = is_fastq_file(input_path)
            .context("Failed to determine if file is FASTQ")?;

        if !is_fasta && !is_fastq {
            return Err(anyhow::anyhow!("Unsupported file format: {}", input_path.display()));
        }

        let file_stem = input_path.file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("output");
            
        let output_name = if format == "fq" {
            format!("{}.fq", file_stem)
        } else {
            format!("{}.fa", file_stem)
        };
        
        let mut output_path = PathBuf::from(sample_output_dir);
        output_path.push(output_name);

        // 根据文件类型处理
        if is_fasta {
            process_fasta(input_path, &output_path, &enzyme, format, input_path.to_string_lossy().ends_with(".gz"))?;
        } else {
            process_fastq(input_path, &output_path, &enzyme, format, input_path.to_string_lossy().ends_with(".gz"))?;
        }
    }

    Ok(())
}



#[allow(dead_code)]
fn is_fasta_file(path: &Path) -> Result<bool> {
    // 检查文件扩展名
    let ext = path.extension()
        .and_then(|s| s.to_str())
        .unwrap_or("")
        .to_lowercase();

    // 如果是压缩文件，获取原始扩展名
    let base_name = path.file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("");
    
    let is_fasta_ext = if base_name.ends_with(".gz") {
        // 移除 .gz 后缀并检查原始扩展名
        let without_gz = base_name.trim_end_matches(".gz");
        without_gz.ends_with(".fa") || 
        without_gz.ends_with(".fasta") || 
        without_gz.ends_with(".fna") || 
        without_gz.ends_with(".ffn") || 
        without_gz.ends_with(".faa") || 
        without_gz.ends_with(".frn")
    } else {
        matches!(ext.as_str(), 
            "fa" | "fasta" | "fna" | "ffn" | "faa" | "frn"
        )
    };

    // 如果扩展名不明确，检查文件内容
    if !is_fasta_ext {
        let mut reader = create_reader(path)?;
        let mut first_char = [0u8; 1];
        if reader.read_exact(&mut first_char).is_ok() {
            return Ok(first_char[0] == b'>');
        }
    }

    Ok(is_fasta_ext)
}

#[allow(dead_code)]
fn is_fastq_file(path: &Path) -> Result<bool> {
    // 检查文件扩展名
    let base_name = path.file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("");
    
    let is_fastq_ext = if base_name.ends_with(".gz") {
        // 移除 .gz 后缀并检查原始扩展名
        let without_gz = base_name.trim_end_matches(".gz");
        without_gz.ends_with(".fq") || 
        without_gz.ends_with(".fastq")
    } else {
        base_name.ends_with(".fq") || 
        base_name.ends_with(".fastq")
    };

    // 如果扩展名不明确，检查文件内容
    if !is_fastq_ext {
        let mut reader = create_reader(path)?;
        let mut first_char = [0u8; 1];
        if reader.read_exact(&mut first_char).is_ok() {
            return Ok(first_char[0] == b'@');
        }
    }

    Ok(is_fastq_ext)
}

// 完全按照sylph方式处理FASTA文件
fn process_fasta_sylph_style(
    input: &Path,
    output: &Path,
    enzyme: &EnzymeSpec,
    format: &str,
    compress: bool,
) -> Result<()> {
    let mut writer = create_writer(output, compress)?;
    let mut stats = ExtractionStats::new();
    
    // 完全按照sylph的模式
    let reader = parse_fastx_file(input);
    if !reader.is_ok() {
        warn!("{} is not a valid fasta/fastq file; skipping.", input.display());
        return Ok(());
    }
    
    let mut reader = reader.unwrap();
    let mut kmer_to_tag_table = FxHashSet::default();
    
    while let Some(record) = reader.next() {
        if record.is_ok() {
            let record = record.expect(&format!("Invalid record for file {} ", input.display()));
            let seq = record.seq();
            let seq_id = String::from_utf8_lossy(record.id());
            
            stats.total_sequences += 1;
            stats.total_sequence_length += seq.len();

            // canonical tag 字节序列（用于写出 FASTA/FASTQ）
            let tags = extract_canonical_tags(&seq, enzyme);

            // 按照sylph的去重模式（现在使用canonical tags）
            for (tag, _len) in tags {
                if kmer_to_tag_table.insert(tag.clone()) {
                    stats.total_tags += 1;
                    write_tags(&mut *writer, &seq_id, &[tag], format)?;
                }
            }
        } else {
            warn!("Invalid record in file {}", input.display());
        }
    }

    log_stats(stats, enzyme);
    Ok(())
}



fn process_fasta(
    input: &Path,
    output: &Path,
    enzyme: &EnzymeSpec,
    format: &str,
    compress: bool,
) -> Result<()> {
    // 直接使用sylph风格的处理
    process_fasta_sylph_style(input, output, enzyme, format, compress)
}

// 按照sylph风格处理FASTQ文件
fn process_fastq_sylph_style(
    input: &Path,
    output: &Path,
    enzyme: &EnzymeSpec,
    format: &str,
    compress: bool,
) -> Result<()> {
    let mut writer = create_writer(output, compress)?;
    let mut stats = ExtractionStats::new();
    
    // 完全按照sylph的模式
    let reader = parse_fastx_file(input);
    if !reader.is_ok() {
        warn!("{} is not a valid fasta/fastq file; skipping.", input.display());
        return Ok(());
    }
    
    let mut reader = reader.unwrap();
    let mut kmer_to_tag_table = FxHashSet::default();
    
    while let Some(record) = reader.next() {
        if record.is_ok() {
            let record = record.expect(&format!("Invalid record for file {} ", input.display()));
            let seq = record.seq();
            let seq_id = String::from_utf8_lossy(record.id());
            
            stats.total_sequences += 1;
            stats.total_sequence_length += seq.len();

            // canonical tag 字节序列（用于写出 FASTA/FASTQ）
            let tags = extract_canonical_tags(&seq, enzyme);

            // 按照sylph的去重模式（现在使用canonical tags）
            for (tag, _len) in tags {
                if kmer_to_tag_table.insert(tag.clone()) {
                    stats.total_tags += 1;
                    write_tags(&mut *writer, &seq_id, &[tag], format)?;
                }
            }
        } else {
            warn!("Invalid record in file {}", input.display());
        }
    }

    log_stats(stats, enzyme);
    Ok(())
}

fn process_fastq(
    input: &Path,
    output: &Path,
    enzyme: &EnzymeSpec,
    format: &str,
    compress: bool,
) -> Result<()> {
    // 直接使用sylph风格的处理
    process_fastq_sylph_style(input, output, enzyme, format, compress)
}

fn write_tags(
    writer: &mut dyn Write,
    seq_id: &str,
    tags: &[TagHash],
    format: &str,
) -> Result<()> {
    for (i, tag) in tags.iter().enumerate() {
        let header = format!("{}_tag{}", seq_id, i + 1);
        match format {
            "fq" => writeln!(writer, "@{}\n{}\n+\n{}", 
                           header, 
                           String::from_utf8_lossy(tag),
                           "~".repeat(tag.len()))
                .context("Failed to write FASTQ record")?,
            _ => writeln!(writer, ">{}\n{}", 
                        header, 
                        String::from_utf8_lossy(tag))
                .context("Failed to write FASTA record")?,
        }
    }
    Ok(())
}

fn create_reader(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path)
        .context(format!("Failed to open input file: {}", path.display()))?;

    // 使用优化的文件大小检测
    let file_size = get_file_size_optimized(path)?;
    let is_compressed = path.to_string_lossy().ends_with(".gz");
    
    // 使用优化的缓冲区大小计算
    let buffer_size = calculate_optimal_buffer_size(file_size, is_compressed);

    Ok(if is_compressed {
        Box::new(BufReader::with_capacity(buffer_size, GzDecoder::new(file)))
    } else {
        Box::new(BufReader::with_capacity(buffer_size, file))
    })
}

fn create_writer(path: &Path, compress: bool) -> Result<Box<dyn Write>> {
    let file = File::create(path)
        .context(format!("Failed to create output file: {}", path.display()))?;

    // 使用优化的缓冲区大小和压缩设置
    let buffer_size = if compress {
        256 * 1024
    } else {
        128 * 1024
    };

    Ok(if compress {
        let compression = get_optimal_compression();
        Box::new(BufWriter::with_capacity(buffer_size, GzEncoder::new(file, compression)))
    } else {
        Box::new(BufWriter::with_capacity(buffer_size, file))
    })
}

#[derive(Debug, Clone)]
pub struct ExtractionStats {
    total_sequences: usize,
    total_tags: usize,
    total_sequence_length: usize,
}



impl ExtractionStats {
    fn new() -> Self {
        Self {
            total_sequences: 0,
            total_tags: 0,
            total_sequence_length: 0,
        }
    }
}

fn log_stats(stats: ExtractionStats, enzyme: &EnzymeSpec) {
    let k = enzyme.patterns[0].as_str().len();
    let total_kmers = if stats.total_sequence_length >= (k - 1) * stats.total_sequences {
        stats.total_sequence_length - (k - 1) * stats.total_sequences
    } else {
        0
    };
    let percentage = calculate_tag_percentage(stats.total_tags, total_kmers);

    // 酶的标签长度在 EnzymeSpec 构造时已查好，直接复用
    let tag_length = enzyme.tag_length;

    let tag_bases_percentage = (stats.total_tags * tag_length) as f64 / stats.total_sequence_length as f64 * 100.0;
    
    println!(
        "\nProcessing complete for {}:\n\
        =============================\n\
        - Total sequences processed: {}\n\
        - Total sequence length: {}\n\
        - Average sequence length: {:.2}\n\
        - Total tags extracted: {}\n\
        - Average tags per sequence: {:.2}\n\
        - Extractable k-mers: {}\n\
        - 2bRAD tag percentage: {:.4}%\n\
        - 2bRAD tag bases percentage: {:.4}%\n\
        - Recognition patterns used: {}",
        enzyme.name,
        stats.total_sequences,
        stats.total_sequence_length,
        stats.total_sequence_length as f32 / stats.total_sequences.max(1) as f32,
        stats.total_tags,
        stats.total_tags as f32 / stats.total_sequences.max(1) as f32,
        total_kmers,
        percentage,
        tag_bases_percentage,
        enzyme.patterns
            .iter()
            .map(|r| r.as_str())
            .collect::<Vec<_>>()
            .join(", ")
    );
}

// 新增函数：处理单对双端测序文件
fn process_paired_fastq_files(
    first_file: &str,
    second_file: &str,
    enzyme: &EnzymeSpec,
    _sample_output_dir: &Path,
    _out_name: Option<&str>,
) -> Result<()> {
    // 从文件名中提取样本名
    let file_stem = Path::new(first_file)
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("unknown")
        .split('.')
        .next()
        .unwrap_or("unknown")
        .to_string();

    // 处理一对文件
    let fa_entries = process_paired_fastq_to_sylsp(
        first_file,
        second_file,
        enzyme,
        &file_stem,
    )?;

    // 注释掉生成单个文件的代码 - 只保留合并后的文件
    // let output_base = Path::new(sample_output_dir).join(&file_stem);
    // let fa_path = output_base.with_extension("fa");
    // let mut fa_writer = create_writer(&fa_path, false)?;

    let mut sylsp_entries = Vec::new();
    for (id, tag, tag_len, sample_source) in &fa_entries {
        let entry = SylspEntry {
            sequence_id: id.clone(),
            tag: *tag,
            tag_length: *tag_len,
            sample_source: sample_source.clone(),
            enzyme: enzyme.name.clone(),
        };
        sylsp_entries.push(entry);
    }

    // 注释掉生成单个sylsp文件的代码
    // let sylsp_path = if let Some(name) = out_name {
    //     Path::new(sample_output_dir).join(format!("{}.sylsp", name))
    // } else {
    //     output_base.with_extension("sylsp")
    // };
    //
    // let sylsp_file = File::create(&sylsp_path)
    //     .context(format!("Failed to create sylsp file: {}", sylsp_path.display()))?;
    // let sylsp_writer = BufWriter::new(sylsp_file);
    // bincode::serialize_into(sylsp_writer, &sylsp_entries)
    //     .context("Failed to serialize sylsp data")?;

    Ok(())
}

pub fn extract(args: ExtractArgs) -> Result<()> {
    // 初始化线程池
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()?;

    // 创建输出目录
    std::fs::create_dir_all(&args.sample_output_dir)
        .context("Failed to create output directory")?;

    // 设置内存限制：未指定时取总物理内存的 75%（旧默认 16GB 会在并行建库时
    // 触发阻塞式 guard 导致死锁，见 check_vram_and_block 注释）
    let max_ram = args.max_ram.unwrap_or_else(default_max_ram_gb);
    if max_ram < 7 {
        return Err(anyhow::anyhow!("Max ram must be >= 7. Exiting."));
    }
    if args.max_ram.is_none() {
        eprintln!("Memory limit: {} GB (default 75% of total RAM; override with --max-ram)", max_ram);
    }

    // 处理单对双端测序文件（-1 和 -2 参数）
    if !args.first_pair.is_empty() && !args.second_pair.is_empty() {
        let enzyme = EnzymeSpec::new(&args.enzyme)?;
        for (first_file, second_file) in args.first_pair.iter().zip(args.second_pair.iter()) {
            safe_process_with_memory_check(max_ram, first_file, || {
                process_paired_fastq_files(
                    first_file,
                    second_file,
                    &enzyme,
                    Path::new(&args.sample_output_dir),
                    args.out_name.as_deref(),
                )
            })?;
        }
    }

    // 处理批处理双端测序文件（--l1 和 --l2 参数）
    if let (Some(first_pair_list), Some(second_pair_list)) = (&args.first_pair_list, &args.second_pair_list) {
        // 读取文件列表
        let first_pairs = read_file_list(first_pair_list)
            .context("Failed to read first pair list")?;
        let second_pairs = read_file_list(second_pair_list)
            .context("Failed to read second pair list")?;

        if first_pairs.len() != second_pairs.len() {
            return Err(anyhow::anyhow!("Number of files in first pair list and second pair list do not match"));
        }

        let enzyme = EnzymeSpec::new(&args.enzyme)?;
        let mut all_sylsp_entries = Vec::new();

        // 并行处理所有配对文件，添加内存监控
        let results: Vec<Result<(String, Vec<SylspEntry>)>> = first_pairs.par_iter()
            .zip(second_pairs.par_iter())
            .map(|(first_file, second_file)| {
                // 检查内存使用
                if let Some(current_memory) = get_memory_usage() {
                    if current_memory > max_ram as f64 {
                        check_vram_and_block(max_ram, first_file);
                    }
                }
                
                let input_path1 = PathBuf::from(first_file);
                let file_stem = input_path1.file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or("unknown")
                    .split('.')
                    .next()
                    .unwrap_or("unknown")
                    .to_string();

                // 处理一对文件
                let fa_entries = process_paired_fastq_to_sylsp(
                    first_file,
                    second_file,
                    &enzyme,
                    &file_stem,
                )?;

                // 注释掉生成单个文件的代码 - 只保留合并后的文件
                // let output_base = Path::new(&args.sample_output_dir).join(&file_stem);
                // let fa_path = output_base.with_extension("fa");
                // let mut fa_writer = create_writer(&fa_path, false)?;

                let mut sylsp_entries = Vec::new();
                for (id, tag, tag_len, sample_source) in &fa_entries {
                    let entry = SylspEntry {
                        sequence_id: id.clone(),
                        tag: *tag,
                        tag_length: *tag_len,
                        sample_source: sample_source.clone(),
                        enzyme: enzyme.name.clone(),
                    };
                    sylsp_entries.push(entry);
                }

                // 注释掉生成单个sylsp文件的代码
                // let sylsp_path = output_base.with_extension("sylsp");
                // let sylsp_file = File::create(&sylsp_path)
                //     .context(format!("Failed to create sylsp file: {}", sylsp_path.display()))?;
                // let sylsp_writer = BufWriter::new(sylsp_file);
                // bincode::serialize_into(sylsp_writer, &sylsp_entries)
                //     .context("Failed to serialize sylsp data")?;

                Ok((file_stem, sylsp_entries))
            })
            .collect();

        // 处理结果并收集所有 sylsp 条目
        for result in results {
            match result {
                Ok((_, entries)) => {
                    all_sylsp_entries.extend(entries);
                },
                Err(e) => eprintln!("Error processing paired files: {}", e),
            }
        }

        // 生成合并的 sylsp 文件
        if !all_sylsp_entries.is_empty() {
            let output_name = args.out_name.as_ref().map_or_else(|| "combined".to_string(), |s| s.clone());
            let combined_sylsp_path = Path::new(&args.sample_output_dir).join(format!("{}.sylsp", output_name));
            let combined_sylsp_file = File::create(&combined_sylsp_path)
                .context(format!("Failed to create combined sylsp file: {}", combined_sylsp_path.display()))?;
            let combined_sylsp_writer = BufWriter::new(combined_sylsp_file);
            
            bincode::serialize_into(combined_sylsp_writer, &all_sylsp_entries)
                .context("Failed to serialize combined sylsp data")?;
        }
    }

    // 处理单端测序文件
    if let Some(read_files) = args.reads {
        // 存储所有 FASTQ 文件的 sylsp 条目
        let mut all_sylsp_entries = Vec::new();
        let enzyme = EnzymeSpec::new(&args.enzyme)?;

        // FASTA 流式写出：边处理边写，避免把所有 tag 字节缓存在内存里。
        let output_name = args.out_name.as_ref().map_or_else(|| "reads".to_string(), |s| s.clone());
        let fa_path = Path::new(&args.sample_output_dir).join(format!("{}.fasta", output_name));
        let mut fa_writer = create_writer(&fa_path, false)?;

        for file in read_files {
            // 检查内存使用
            if let Some(current_memory) = get_memory_usage() {
                if current_memory > max_ram as f64 {
                    check_vram_and_block(max_ram, &file);
                }
            }

            let input_path = PathBuf::from(&file);
            let file_stem = input_path.file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("unknown")
                .split('.')
                .next()
                .unwrap_or("unknown")
                .to_string();

            let reader = fastq::Reader::new(create_reader(&input_path)?);
            let mut stats = ExtractionStats::new();

            for result in reader.records() {
                let record = result.context("Failed to read FASTQ record")?;
                stats.total_sequences += 1;
                stats.total_sequence_length += record.seq().len();

                let tags = extract_canonical_tags(record.seq(), &enzyme);

                for (i, (tag, tag_len)) in tags.iter().enumerate() {
                    let id = format!("{}_tag{}", record.id(), i + 1);
                    writeln!(fa_writer, ">{}\n{}", id, String::from_utf8_lossy(tag))
                        .context("Failed to write FASTA record")?;

                    all_sylsp_entries.push(SylspEntry {
                        sequence_id: id,
                        tag: hash_bytes(tag),
                        tag_length: *tag_len,
                        sample_source: file_stem.clone(),
                        enzyme: enzyme.name.clone(),
                    });
                }

                stats.total_tags += tags.len();
            }

            log_stats(stats, &enzyme);
        }

        // 生成 .sylsp 文件
        let sylsp_path = Path::new(&args.sample_output_dir).join(format!("{}.sylsp", output_name));
        let sylsp_file = File::create(&sylsp_path)
            .context(format!("Failed to create sylsp file: {}", sylsp_path.display()))?;
        let sylsp_writer = BufWriter::new(sylsp_file);

        bincode::serialize_into(sylsp_writer, &all_sylsp_entries)
            .context("Failed to serialize sylsp data")?;
    }

    // 处理基因组列表文件
    if let Some(genome_list) = &args.genome_list {
        let file = File::open(genome_list)
            .context(format!("Failed to open genome list file: {}", genome_list))?;
        let reader = BufReader::new(file);
        let genome_files: Vec<String> = reader.lines()
            .filter_map(|line| line.ok())
            .collect();

        let enzyme = EnzymeSpec::new(&args.enzyme)?;
        let mut all_syldb_entries = Vec::new();
        
        // 并行处理所有 FASTA 文件，添加内存监控
        let results: Vec<Result<Vec<SyldbEntry>>> = genome_files.par_iter()
            .map(|file| {
                // 检查内存使用
                if let Some(current_memory) = get_memory_usage() {
                    if current_memory > max_ram as f64 {
                        check_vram_and_block(max_ram, file);
                    }
                }
                
                let input_path = Path::new(file);
                let output_base = Path::new(&args.sample_output_dir).join(input_path.file_stem().unwrap_or_default());
                process_fasta_to_syldb(
                    input_path,
                    &output_base,
                    &enzyme,
                    &args.format,
                    file.ends_with(".gz"),
                    !args.no_tag_seqs,
                )
            })
            .collect();

        // 收集所有结果
        for (file, result) in genome_files.iter().zip(results) {
            match result {
                Ok(mut entries) => {
                    // 为每个条目添加基因组来源信息
                    for entry in &mut entries {
                        entry.genome_source = file.clone();
                    }
                    all_syldb_entries.extend(entries);
                },
                Err(e) => {
                    eprintln!("Error processing FASTA file: {}", e);
                }
            }
        }

        // 生成合并的 .syldb 文件
        if !all_syldb_entries.is_empty() {
            let output_name = args.out_name.as_ref().map_or_else(|| "combined".to_string(), |s| s.clone());
            let combined_syldb_path = Path::new(&args.sample_output_dir).join(format!("{}.syldb", output_name));
            let combined_syldb_file = File::create(&combined_syldb_path)
                .context(format!("Failed to create combined syldb file: {}", combined_syldb_path.display()))?;
            let combined_syldb_writer = BufWriter::new(combined_syldb_file);
            
            bincode::serialize_into(combined_syldb_writer, &all_syldb_entries)
                .context("Failed to serialize combined syldb data")?;
        }
    }

    // 处理基因组文件
    if let Some(genome_files) = &args.genomes {
        let enzyme = EnzymeSpec::new(&args.enzyme)?;
        let mut all_syldb_entries = Vec::new();
        
        // 并行处理所有 FASTA 文件，添加内存监控
        let results: Vec<Result<Vec<SyldbEntry>>> = genome_files.par_iter()
            .map(|file| {
                // 检查内存使用
                if let Some(current_memory) = get_memory_usage() {
                    if current_memory > max_ram as f64 {
                        check_vram_and_block(max_ram, file);
                    }
                }
                
                let input_path = Path::new(file);
                let output_base = Path::new(&args.sample_output_dir).join(input_path.file_stem().unwrap_or_default());
                process_fasta_to_syldb(
                    input_path,
                    &output_base,
                    &enzyme,
                    &args.format,
                    file.ends_with(".gz"),
                    !args.no_tag_seqs,
                )
            })
            .collect();

        // 收集所有结果
        for (file, result) in genome_files.iter().zip(results) {
            match result {
                Ok(mut entries) => {
                    // 为每个条目添加基因组来源信息
                    for entry in &mut entries {
                        entry.genome_source = file.clone();
                    }
                    all_syldb_entries.extend(entries);
                },
                Err(e) => {
                    eprintln!("Error processing FASTA file: {}", e);
                }
            }
        }

        // 生成合并的 .syldb 文件
        if !all_syldb_entries.is_empty() {
            let output_name = args.out_name.as_ref().map_or_else(|| "combined".to_string(), |s| s.clone());
            let combined_syldb_path = Path::new(&args.sample_output_dir).join(format!("{}.syldb", output_name));
            let combined_syldb_file = File::create(&combined_syldb_path)
                .context(format!("Failed to create combined syldb file: {}", combined_syldb_path.display()))?;
            let combined_syldb_writer = BufWriter::new(combined_syldb_file);
            
            bincode::serialize_into(combined_syldb_writer, &all_syldb_entries)
                .context("Failed to serialize combined syldb data")?;
        }
    }

    // 处理样本列表文件
    if let Some(sample_list) = &args.sample_list {
        let mut all_sylsp_entries = Vec::new();
        let enzyme = EnzymeSpec::new(&args.enzyme)?;
        
        // 读取样本列表文件
        let file = File::open(sample_list)
            .context(format!("Failed to open sample list file: {}", sample_list))?;
        let reader = BufReader::new(file);
        
        // 并行处理所有样本文件
        let sample_files: Vec<String> = reader.lines()
            .filter_map(|line| line.ok())
            .collect();
            
        // 使用FxHashMap优化样本处理
        let sample_stats = Arc::new(Mutex::new(SampleStatsMap::default()));
        
        let results: Vec<Result<(String, Vec<SylspEntry>)>> = sample_files.par_iter()
            .map(|file| {
                // 检查内存使用
                if let Some(current_memory) = get_memory_usage() {
                    if current_memory > max_ram as f64 {
                        check_vram_and_block(max_ram, file);
                    }
                }
                
                let input_path = PathBuf::from(file);
                // 修正样本名提取逻辑
                let file_name = input_path.file_name()
                    .and_then(|s| s.to_str())
                    .unwrap_or("unknown");
                let file_stem = file_name.split('.').next().unwrap_or("unknown").to_string();
                
                let reader = fastq::Reader::new(create_reader(&input_path)?);
                let mut sylsp_entries = Vec::new();
                let mut stats = ExtractionStats::new();

                for result in reader.records() {
                    let record = result.context("Failed to read FASTQ record")?;
                    stats.total_sequences += 1;
                    stats.total_sequence_length += record.seq().len();

                    let tags = extract_tag_hashes(record.seq(), &enzyme);
                    stats.total_tags += tags.len();

                    for (i, (tag, tag_len)) in tags.iter().enumerate() {
                        sylsp_entries.push(SylspEntry {
                            sequence_id: format!("{}_tag{}", record.id(), i + 1),
                            tag: *tag,
                            tag_length: *tag_len,
                            sample_source: file_stem.clone(), // 用文件名去除扩展名作为样本名
                            enzyme: enzyme.name.clone(),
                        });
                    }
                }

                // 更新全局统计
                let mut global_stats = sample_stats.lock().unwrap();
                global_stats.insert(file_stem.clone(), stats.clone());

                log_stats(stats, &enzyme);
                Ok((file_stem, sylsp_entries))
            })
            .collect();
            
        // 处理每个样本的结果
        for result in results {
            match result {
                Ok((_file_stem, sylsp_entries)) => {
                    // 注释掉为每个样本生成独立文件的代码
                    // let fa_path = Path::new(&args.sample_output_dir)
                    //     .join(format!("{}.fasta", file_stem));
                    // let mut fa_writer = create_writer(&fa_path, false)?;
                    // 
                    // for (id, tag) in fa_entries {
                    //     writeln!(fa_writer, ">{}\n{}", id, String::from_utf8_lossy(&tag))
                    //         .context("Failed to write FASTA record")?;
                    // }
                    // 
                    // // 为每个样本生成独立的 sylsp 文件
                    // let sample_sylsp_path = Path::new(&args.sample_output_dir)
                    //     .join(format!("{}.sylsp", file_stem));
                    // let sample_sylsp_file = File::create(&sample_sylsp_path)
                    //     .context(format!("Failed to create sylsp file: {}", sample_sylsp_path.display()))?;
                    // let sample_sylsp_writer = BufWriter::new(sample_sylsp_file);
                    // 
                    // bincode::serialize_into(sample_sylsp_writer, &sylsp_entries)
                    //     .context(format!("Failed to serialize sylsp data for sample: {}", file_stem))?;
                    
                    // 收集所有 sylsp 条目用于合并
                    all_sylsp_entries.extend(sylsp_entries);
                },
                Err(e) => eprintln!("Error processing sample file: {}", e),
            }
        }
        
        // 生成合并的 .sylsp 文件
        let output_name = args.out_name.as_ref().map_or_else(|| "combined".to_string(), |s| s.clone());
        let sylsp_path = Path::new(&args.sample_output_dir).join(format!("{}.sylsp", output_name));
        let sylsp_file = File::create(&sylsp_path)
            .context(format!("Failed to create combined sylsp file: {}", sylsp_path.display()))?;
        let sylsp_writer = BufWriter::new(sylsp_file);
        
        bincode::serialize_into(sylsp_writer, &all_sylsp_entries)
            .context("Failed to serialize combined sylsp data")?;
    }

    Ok(())
}

fn process_fasta_to_syldb(
    input: &Path,
    _output_base: &Path,
    enzyme: &EnzymeSpec,
    _format: &str,
    _compress: bool,
    store_tag_seqs: bool,
) -> Result<Vec<SyldbEntry>> {
    let enzyme_name = enzyme.name.clone();
    // 注释掉生成单个.fa文件的代码
    // let fa_path = output_base.with_extension("fa");
    // let mut fa_writer = BufWriter::with_capacity(64 * 1024, File::create(&fa_path)?);
    
    let mut stats = ExtractionStats::new();
    // 预分配容量 - 估计每个序列平均产生50个标签
    let mut syldb_entries = Vec::with_capacity(100);

    // 读取和处理 FASTA 记录
    let reader = create_reader(input)?;
    for record in fasta::Reader::new(reader).records() {
        let record = record.context("Failed to read FASTA record")?;
        let seq_len = record.seq().len();
        stats.total_sequences += 1;
        stats.total_sequence_length += seq_len;

        // 提取 canonical tag 字节序列及其哈希；保留序列以支持 error-tolerant matching。
        let tag_items = extract_canonical_tags(record.seq(), enzyme);
        stats.total_tags += tag_items.len();
        let mut tags = Vec::with_capacity(tag_items.len());
        let mut tag_lengths = Vec::with_capacity(tag_items.len());
        let mut tag_seqs = if store_tag_seqs {
            Some(Vec::with_capacity(tag_items.len()))
        } else {
            None
        };
        for (tag, tag_len) in tag_items {
            tags.push(hash_bytes(&tag));
            tag_lengths.push(tag_len);
            if let Some(seqs) = tag_seqs.as_mut() {
                seqs.push(tag);
            }
        }

        // 创建 syldb 条目
        let entry = SyldbEntry {
            sequence_id: record.id().to_string(),
            tags,
            tag_lengths,
            genome_source: input.to_string_lossy().to_string(),
            tag_uniqueness: None, // 初始时未标记，将由mark命令处理
            tag_seqs,
            enzyme: enzyme_name.clone(),
        };
        syldb_entries.push(entry);
    }

    // 注释掉生成单个.syldb文件的代码
    // let syldb_path = output_base.with_extension("syldb");
    // let syldb_file = File::create(&syldb_path)
    //     .context(format!("Failed to create syldb file: {}", syldb_path.display()))?;
    // let syldb_writer = BufWriter::with_capacity(64 * 1024, syldb_file);
    // 
    // // 使用标准序列化 API
    // bincode::serialize_into(syldb_writer, &syldb_entries)
    //     .context("Failed to serialize syldb data")?;

    log_stats(stats, enzyme);
    Ok(syldb_entries)
}



fn process_paired_fastq_to_sylsp(
    input1: &str,
    input2: &str,
    enzyme: &EnzymeSpec,
    sample_source: &str,
) -> Result<Vec<(String, Hash, u8, String)>> {
    let reader1 = fastq::Reader::new(create_reader(Path::new(input1))?);
    let reader2 = fastq::Reader::new(create_reader(Path::new(input2))?);
    let mut stats = ExtractionStats::new();
    let mut entries = Vec::new();

    let mut iter1 = reader1.records();
    let mut iter2 = reader2.records();

    loop {
        let record1 = match iter1.next() {
            Some(Ok(r)) => r,
            Some(Err(e)) => return Err(anyhow::anyhow!("Error reading first pair: {}", e)),
            None => break,
        };

        let record2 = match iter2.next() {
            Some(Ok(r)) => r,
            Some(Err(e)) => return Err(anyhow::anyhow!("Error reading second pair: {}", e)),
            None => break,
        };

        let seq_len1 = record1.seq().len();
        let seq_len2 = record2.seq().len();
        stats.total_sequences += 1;
        stats.total_sequence_length += seq_len1 + seq_len2;

        // 每条 read 内 tag 已由 extract_tag_hashes 去重；read id + tag index 天然唯一，
        // 原先的 seen_pairs 集合永远不会命中重复，纯属浪费，已移除。
        let tags1 = extract_tag_hashes(record1.seq(), enzyme);
        let tags2 = extract_tag_hashes(record2.seq(), enzyme);
        stats.total_tags += tags1.len() + tags2.len();

        for (i, (tag, tag_len)) in tags1.iter().enumerate() {
            entries.push((format!("{}_{}", record1.id(), i + 1), *tag, *tag_len, sample_source.to_string()));
        }
        for (i, (tag, tag_len)) in tags2.iter().enumerate() {
            entries.push((format!("{}_{}", record2.id(), i + 1), *tag, *tag_len, sample_source.to_string()));
        }
    }

    log_stats(stats, enzyme);
    Ok(entries)
}

fn calculate_tag_percentage(tag_count: usize, total_kmers: usize) -> f64 {
    if total_kmers == 0 {
        0.0
    } else {
        (tag_count as f64 / total_kmers as f64) * 100.0
    }
}



fn read_file_list(path: &str) -> Result<Vec<String>> {
    let file = File::open(path)
        .context(format!("Failed to open file list: {}", path))?;
    let reader = BufReader::new(file);
    let files: Vec<String> = reader.lines()
        .filter_map(|line| line.ok())
        .collect();
    Ok(files)
}




#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_process_fasta_to_syldb_tag_seqs_optional() {
        let dir = std::env::temp_dir();
        let fa_path = dir.join("sylph_test_tag_seqs_optional.fa");
        std::fs::write(&fa_path, ">seq1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n").unwrap();
        let out_base = dir.join("sylph_test_tag_seqs_optional");
        let enzyme = EnzymeSpec::new("BcgI").unwrap();

        let with_seqs = process_fasta_to_syldb(&fa_path, &out_base, &enzyme, "fa", false, true).unwrap();
        assert_eq!(with_seqs.len(), 1);
        assert!(with_seqs[0].tag_seqs.is_some());

        let without_seqs = process_fasta_to_syldb(&fa_path, &out_base, &enzyme, "fa", false, false).unwrap();
        assert_eq!(without_seqs.len(), 1);
        assert!(without_seqs[0].tag_seqs.is_none());
        // 哈希与长度不受 flag 影响
        assert_eq!(with_seqs[0].tags, without_seqs[0].tags);
        assert_eq!(with_seqs[0].tag_lengths, without_seqs[0].tag_lengths);

        std::fs::remove_file(&fa_path).ok();
    }

    #[test]
    fn test_one_mismatch_canonical_hashes() {
        // 简单 tag：AAA -> canonical 仍是 AAA（因为 TTT revcomp > AAA）
        let tag = b"AAA";
        let hashes = one_mismatch_canonical_hashes(tag);
        // 3 positions * 3 alts = 9 neighbors
        assert_eq!(hashes.len(), 9);

        // 所有 neighbor hash 都应该与 exact hash 不同
        let exact = hash_bytes(tag);
        assert!(!hashes.contains(&exact));

        // 检查一个特定 neighbor：AAC 的 canonical 是 AAC（GGTTT revcomp > AAC）
        let expected = hash_bytes(b"AAC");
        assert!(hashes.contains(&expected));
    }

    #[test]
    fn test_canonical_neighbor_consistency() {
        // 样本中提取到的 tag 是 canonical 形式；reference 的 1-mismatch neighbor 也必须是 canonical 形式才能匹配。
        // 例：reference tag = "AAACCC"，样本 tag = "AAACCG"（1 mismatch）。
        // canonical("AAACCG") = "AAACCG"（因为 revcomp = "CGGGTT" > "AAACCG"）。
        let ref_tag = b"AAACCC";
        let sample_tag = b"AAACCG";
        let neighbors = one_mismatch_canonical_hashes(ref_tag);
        let sample_hash = hash_bytes(sample_tag);
        assert!(neighbors.contains(&sample_hash));
    }

    #[test]
    fn test_p_detect_one_mismatch() {
        // a=1.0 时必然检出
        assert!((p_detect_one_mismatch(1.0, 30) - 1.0).abs() < 1e-12);
        // a=0.0 时只能由 mismatch 检出，但 1 mismatch 也需要至少 ℓ-1 个匹配，所以 a=0 时仍是 0
        assert!((p_detect_one_mismatch(0.0, 30) - 0.0).abs() < 1e-12);
        // a=0.95, len=30：应略高于 exact
        let p_exact = p_detect_exact(0.95, 30);
        let p_mm = p_detect_one_mismatch(0.95, 30);
        assert!(p_mm > p_exact);
    }

    #[test]
    fn test_ani_from_containment_one_mismatch() {
        // containment = 1.0 -> ANI = 1.0
        assert!((ani_from_containment_one_mismatch(1.0, 30) - 1.0).abs() < 1e-12);
        // containment = 0.0 -> ANI = 0.0
        assert!((ani_from_containment_one_mismatch(0.0, 30) - 0.0).abs() < 1e-12);
        // 允许 1 mismatch 时，同样的 containment 对应更低的真实 ANI（因为检出概率更高）
        let c: f64 = 0.8;
        let a_exact = c.powf(1.0 / 30.0);
        let a_mm = ani_from_containment_one_mismatch(c, 30);
        assert!(a_mm < a_exact, "a_mm={} should be < a_exact={}", a_mm, a_exact);
        // 自检：p_detect_one_mismatch(a_mm, 30) 应接近 c
        let recovered = p_detect_one_mismatch(a_mm, 30);
        assert!((recovered - c).abs() < 1e-6, "recovered={} != c={}", recovered, c);
    }
}
