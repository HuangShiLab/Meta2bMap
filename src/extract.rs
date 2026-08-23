use anyhow::{Context, Result};
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
    /// 单趟扫描引擎（按 pattern 数量分派，见 EnzymeSpec::new 的注释）。
    pub scanner: ScanMode,
    /// 每个 pattern 对应的 tag 长度。多酶模式下 pattern 来自不同酶，长度可能不同。
    pub pattern_tag_lengths: Vec<usize>,
    /// 主 tag 长度（日志/统计使用，取第一个酶的 tag 长度）。
    pub tag_length: usize,
}

/// 单趟扫描引擎。两种实现产出的命中集合与顺序完全一致（均由单元测试对照
/// 旧版逐 pattern 扫描验证），仅性能特征不同：
/// - Regex：组合 alternation 共享各分支的字面量 prefilter，pattern 少时 lazy DFA
///   装得下，单趟最快；pattern 多时 DFA 状态爆炸退化为 PikeVM 逐字节模拟，极慢。
/// - Ac：Aho-Corasick 对所有 pattern 的最长字面量核心串做 overlapping 扫描得到
///   候选起点，再做逐位置掩码验证。overlapping 迭代没有 prefilter，pattern 少时
///   不如 Regex 路径，但随 pattern 数扩展良好（`--enzyme all` 29 个 pattern）。
#[derive(Debug)]
pub enum ScanMode {
    Regex {
        /// 所有 pattern 的组合 alternation，单趟定位任一 pattern 的下一个命中起点。
        combined: Regex,
        /// 每个 pattern 的锚定版本（^(?:pat)），在命中起点确认具体哪些 pattern 匹配。
        anchored: Vec<Regex>,
    },
    Ac {
        /// 所有核心串组成的 Aho-Corasick 自动机（pattern id 与 `patterns` 索引一一对应，
        /// 重复核心串各自占位）。
        ac: aho_corasick::AhoCorasick,
        /// 每个 pattern 的逐位置碱基掩码（bit: A=1,C=2,G=4,T=8），
        /// 用于命中窗口的直接验证（等价于原正则的定长匹配）。
        masks: Vec<Vec<u8>>,
        /// 每个 pattern 选定的触发核心串（最长字面量段）在 pattern 内的偏移。
        core_off: Vec<usize>,
        /// 每个 pattern 的第二长字面量段（偏移, 字节），全掩码验证前的快速过滤。
        core2: Vec<Option<(usize, Vec<u8>)>>,
    },
}

/// 少于等于该数量的 pattern 走 Regex 路径（实测 BcgI 2 个 pattern 时比 AC 快约 2x）；
/// 超过则走 AC（29 个 pattern 时 Regex 路径比基线慢几十倍）。
const REGEX_SCAN_MAX_PATTERNS: usize = 8;

/// 碱基到掩码位的查表（A=1,C=2,G=4,T=8，其他为 0）。
const BASE_BIT: [u8; 256] = {
    let mut table = [0u8; 256];
    table[b'A' as usize] = 1;
    table[b'C' as usize] = 2;
    table[b'G' as usize] = 4;
    table[b'T' as usize] = 8;
    table
};

/// 把酶切正则（仅由 ACGT 字面量、[...] 字符类和定长 {n} 重复组成）解析为
/// 逐位置掩码 + 字面量段列表（偏移, 字节）。语法超出的 pattern 直接报错。
fn parse_pattern_masks(pat: &str) -> Result<(Vec<u8>, Vec<(usize, Vec<u8>)>)> {
    let b = pat.as_bytes();
    let mut i = 0usize;
    let mut masks = Vec::new();
    let mut runs: Vec<(usize, Vec<u8>)> = Vec::new();
    while i < b.len() {
        match b[i] {
            b'[' => {
                let j = b[i..]
                    .iter()
                    .position(|&c| c == b']')
                    .map(|p| p + i)
                    .ok_or_else(|| anyhow::anyhow!("Unclosed class in pattern: {}", pat))?;
                let mut m = 0u8;
                for &c in &b[i + 1..j] {
                    let bit = BASE_BIT[c as usize];
                    if bit == 0 {
                        return Err(anyhow::anyhow!("Unsupported char in pattern: {}", pat));
                    }
                    m |= bit;
                }
                i = j + 1;
                // 可选的 {n} 定长重复
                let rep = if i < b.len() && b[i] == b'{' {
                    let k = b[i..]
                        .iter()
                        .position(|&c| c == b'}')
                        .map(|p| p + i)
                        .ok_or_else(|| anyhow::anyhow!("Unclosed quantifier in pattern: {}", pat))?;
                    let n: usize = std::str::from_utf8(&b[i + 1..k])
                        .ok()
                        .and_then(|s| s.parse().ok())
                        .ok_or_else(|| anyhow::anyhow!("Bad quantifier in pattern: {}", pat))?;
                    i = k + 1;
                    n
                } else {
                    1
                };
                for _ in 0..rep {
                    masks.push(m);
                }
            }
            c @ (b'A' | b'C' | b'G' | b'T') => {
                let _ = c;
                let start = i;
                let off = masks.len();
                while i < b.len() && matches!(b[i], b'A' | b'C' | b'G' | b'T') {
                    masks.push(BASE_BIT[b[i] as usize]);
                    i += 1;
                }
                if i < b.len() && b[i] == b'{' {
                    // 当前酶表不存在「字面量段后接 {n}」的情况，拒绝以避免静默错解
                    return Err(anyhow::anyhow!(
                        "Quantifier on literal run unsupported: {}",
                        pat
                    ));
                }
                runs.push((off, b[start..i].to_vec()));
            }
            _ => return Err(anyhow::anyhow!("Unsupported syntax in pattern: {}", pat)),
        }
    }
    Ok((masks, runs))
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
        let mut all_pattern_strs = Vec::new();
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
                all_pattern_strs.push(*pat);
                all_lengths.push(tag_length);
            }
        }

        // 单趟扫描引擎按 pattern 数量分派（见 ScanMode 注释）：
        // 少量 pattern 用组合 alternation（共享 prefilter，实测最快）；
        // 大量 pattern 用 AC 核心串触发 + 掩码验证（随 pattern 数扩展良好）。
        // 曾尝试对 29 个 pattern 也用 alternation + 放大 DFA 缓存，仍会退化，比
        // 逐 pattern 基线慢几十倍。
        let scanner = if all_pattern_strs.len() <= REGEX_SCAN_MAX_PATTERNS {
            let combined = regex::bytes::RegexBuilder::new(&format!(
                "(?:{})",
                all_pattern_strs.join("|")
            ))
            .dfa_size_limit(64 * (1 << 20))
            .build()
            .context("Invalid combined enzyme pattern")?;
            let mut anchored = Vec::with_capacity(all_pattern_strs.len());
            for pat in &all_pattern_strs {
                anchored.push(
                    Regex::new(&format!("^(?:{})", pat))
                        .context(format!("Invalid anchored regex pattern: {}", pat))?,
                );
            }
            ScanMode::Regex { combined, anchored }
        } else {
            let mut all_masks = Vec::with_capacity(all_pattern_strs.len());
            let mut core_off = Vec::with_capacity(all_pattern_strs.len());
            let mut core2: Vec<Option<(usize, Vec<u8>)>> =
                Vec::with_capacity(all_pattern_strs.len());
            let mut run_store: Vec<Vec<(usize, Vec<u8>)>> =
                Vec::with_capacity(all_pattern_strs.len());
            for pat in &all_pattern_strs {
                let (masks, mut runs) = parse_pattern_masks(pat)?;
                if runs.is_empty() {
                    return Err(anyhow::anyhow!("Pattern has no literal core: {}", pat));
                }
                // 按长度降序：最长段为触发核心，次长段（>=2bp）为快速过滤器
                runs.sort_by_key(|(_, s)| std::cmp::Reverse(s.len()));
                let c2 = if runs.len() > 1 && runs[1].1.len() >= 2 {
                    Some(runs[1].clone())
                } else {
                    None
                };
                core_off.push(runs[0].0);
                core2.push(c2);
                all_masks.push(masks);
                run_store.push(runs);
            }
            let core_strs: Vec<&[u8]> = run_store.iter().map(|r| r[0].1.as_slice()).collect();
            let ac = aho_corasick::AhoCorasick::builder()
                .match_kind(aho_corasick::MatchKind::Standard)
                .kind(Some(aho_corasick::AhoCorasickKind::DFA))
                .build(core_strs)
                .map_err(|e| anyhow::anyhow!("Failed to build AC automaton: {}", e))?;
            ScanMode::Ac {
                ac,
                masks: all_masks,
                core_off,
                core2,
            }
        };

        let primary_length = all_lengths[0];
        let display_name = names.join(",");

        Ok(Self {
            name: display_name,
            patterns: all_patterns,
            scanner,
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

/// 每条 read/contig 复用的临时缓冲区，在文件的 record 循环外创建一次、
/// 循环内 clear 复用，避免 per-read 的 Vec/FxHashSet 堆分配。
struct TagBufs {
    /// find_all_tag_positions_into 的原始命中 (pattern_idx, start, len)，排序前的暂存。
    hits: Vec<(usize, usize, usize)>,
    /// 排序后的 (start, len)，顺序与旧版逐 pattern 扫描完全一致。
    positions: Vec<(usize, usize)>,
    /// 按 u64 哈希去重。
    seen: FxHashSet<Hash>,
    /// canonical 化的栈上缓冲。
    cbuf: [u8; MAX_TAG_LEN],
}

impl Default for TagBufs {
    fn default() -> Self {
        Self {
            hits: Vec::new(),
            positions: Vec::new(),
            seen: FxHashSet::default(),
            cbuf: [0u8; MAX_TAG_LEN],
        }
    }
}

/// `extract_tag_hashes` 的缓冲复用版本：结果写入 `out`（内部先 clear）。
fn extract_tag_hashes_into(
    seq: &[u8],
    enzyme: &EnzymeSpec,
    bufs: &mut TagBufs,
    out: &mut Vec<(Hash, u8)>,
) {
    out.clear();
    bufs.seen.clear();
    find_all_tag_positions_into(seq, enzyme, &mut bufs.hits, &mut bufs.positions);
    for &(offset, len) in &bufs.positions {
        let n = canonicalize_into(&seq[offset..offset + len], &mut bufs.cbuf);
        let h = hash_bytes(&bufs.cbuf[..n]);
        if bufs.seen.insert(h) {
            out.push((h, n as u8));
        }
    }
}

/// `extract_canonical_tags` 的缓冲复用版本，同时返回去重用的哈希，
/// 避免调用方对同一 tag 二次 hash_bytes。结果写入 `out`（内部先 clear）。
fn extract_canonical_tags_into(
    seq: &[u8],
    enzyme: &EnzymeSpec,
    bufs: &mut TagBufs,
    out: &mut Vec<(Hash, TagHash, u8)>,
) {
    out.clear();
    bufs.seen.clear();
    find_all_tag_positions_into(seq, enzyme, &mut bufs.hits, &mut bufs.positions);
    for &(offset, len) in &bufs.positions {
        let n = canonicalize_into(&seq[offset..offset + len], &mut bufs.cbuf);
        let h = hash_bytes(&bufs.cbuf[..n]);
        if bufs.seen.insert(h) {
            out.push((h, bufs.cbuf[..n].to_vec(), n as u8));
        }
    }
}

/// `extract_canonical_tags_into` 的位置感知版本，额外返回每个（去重后）tag 在序列上的
/// bp 偏移（首次出现位置）。仅数据库构建路径使用，样本路径不需要位置信息。
fn extract_canonical_tags_pos_into(
    seq: &[u8],
    enzyme: &EnzymeSpec,
    bufs: &mut TagBufs,
    out: &mut Vec<(Hash, TagHash, u8, u32)>,
) {
    out.clear();
    bufs.seen.clear();
    find_all_tag_positions_into(seq, enzyme, &mut bufs.hits, &mut bufs.positions);
    for &(offset, len) in &bufs.positions {
        let n = canonicalize_into(&seq[offset..offset + len], &mut bufs.cbuf);
        let h = hash_bytes(&bufs.cbuf[..n]);
        if bufs.seen.insert(h) {
            out.push((h, bufs.cbuf[..n].to_vec(), n as u8, offset as u32));
        }
    }
}

/// rust-bio fastq 的 `id()` 只在第一个空格处截断（tab 不截断，见 bio-1.6 fastq.rs
/// `splitn(2, ' ')`），而 needletail 的 `id()` 返回完整 header 行。
/// 统一按首个空格截断，保证切换解析器后 id 逐字节一致。
#[inline]
fn fastq_id(id: &[u8]) -> &[u8] {
    match id.iter().position(|&b| b == b' ') {
        Some(p) => &id[..p],
        None => id,
    }
}

/// rust-bio fasta 的 `id()` 在第一个任意空白字符处截断
/// （`splitn(2, char::is_whitespace)`），needletail 返回完整 header 行。
/// 统一按首个 ASCII 空白截断，与基线一致。
#[inline]
fn fasta_id(id: &[u8]) -> &[u8] {
    match id.iter().position(|b| b.is_ascii_whitespace()) {
        Some(p) => &id[..p],
        None => id,
    }
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
/// This reproduces the same hit set and ordering in a *single* combined
/// scan (see `find_all_tag_positions_into`, which dispatches between two
/// equivalent single-pass engines based on pattern count). All enzyme
/// patterns use fixed-count `{n}` repetitions only (no variable-length
/// quantifiers), so every match is guaranteed to have length
/// `enzyme.tag_length` regardless of where it starts — no separate
/// anchoring or length check is needed for the hit length itself.
#[cfg(test)]
fn find_all_tag_positions(seq: &[u8], enzyme: &EnzymeSpec) -> Vec<(usize, usize)> {
    let mut hits = Vec::new();
    let mut out = Vec::new();
    find_all_tag_positions_into(seq, enzyme, &mut hits, &mut out);
    out
}

/// 单趟扫描版本，按 `enzyme.scanner` 分派两种等价实现：
///
/// - Regex：组合 alternation 的 `find_at` 定位「任一 pattern 的下一个命中起点」，
///   然后在该起点用各 pattern 的锚定正则确认具体哪些 pattern 命中。
/// - Ac：Aho-Corasick 对所有 pattern 的最长字面量核心串做一遍 overlapping 扫描
///   得到候选起点，再用逐位置碱基掩码做定长窗口验证。
///
/// 与旧版「逐 pattern 全序列扫描 + 命中后 rewind 到 match_start+1」的等价性论证：
/// - 旧版对 pattern i 报告所有满足「pattern i 在起点 s 处完整匹配」的 s（正则定长，
///   从 s+1 继续找不会遗漏任何起点），即命中集合为 {(s, i): pattern i 匹配于 s}。
/// - Regex 路径：游标每次只前进到「上一个命中起点+1」。若某起点 s 处有任一 pattern
///   命中且尚未被访问，则 `find_at(seq, cursor)`（cursor ≤ s）返回的最左命中起点 s'
///   满足 cursor ≤ s' ≤ s；归纳可知所有命中起点都会按升序被访问到，唯一被「遮蔽」的
///   情况是同一起点多个 pattern 同时命中——这正是锚定逐一确认要补回的部分。
/// - Ac 路径：pattern i 在 s 处完整匹配 ⟹ 其最长字面量核心串必在 s+off_i 处出现
///   ⟹ AC overlapping 扫描必然报告该核心命中（overlapping 迭代不遗漏任何位置）。
///   候选验证用的是与原正则逐位置完全等价的碱基掩码（[ACGT]/[CT] 等字符类展开为
///   位掩码，字面量为单 bit），通过验证 ⟺ 原正则在 s 处完整匹配。
///
/// 因此两条路径得到的命中集合都与旧版完全一致；同一起点多个 pattern 同时命中的
/// 情况天然支持（每个 pattern 各自触发、各自验证）。
///
/// 输出顺序：旧版是「按 pattern 分组、组内按起点升序」。两条路径的原始命中都按
/// 起点升序（AC 按终点升序，但同一 pattern 核心串定长，组内等价于起点升序）、
/// pattern 交错；按 (pattern_idx, start) 排序后与旧版逐字节一致（下游按哈希去重
/// 保留首次出现顺序，顺序必须保持一致）。
fn find_all_tag_positions_into(
    seq: &[u8],
    enzyme: &EnzymeSpec,
    hits: &mut Vec<(usize, usize, usize)>,
    out: &mut Vec<(usize, usize)>,
) {
    hits.clear();
    out.clear();
    match &enzyme.scanner {
        ScanMode::Regex { combined, anchored } => {
            let mut pos = 0usize;
            while pos <= seq.len() {
                match combined.find_at(seq, pos) {
                    Some(m) => {
                        let mstart = m.start();
                        let tail = &seq[mstart..];
                        for (idx, (anch, &tag_len)) in
                            anchored.iter().zip(&enzyme.pattern_tag_lengths).enumerate()
                        {
                            if anch.is_match(tail) {
                                hits.push((idx, mstart, tag_len));
                            }
                        }
                        pos = mstart + 1; // rewind: mirrors Perl's `pos($seq) = match_start + 1`
                    }
                    None => break,
                }
            }
        }
        ScanMode::Ac {
            ac,
            masks,
            core_off,
            core2,
        } => {
            for m in ac.find_overlapping_iter(seq) {
                let idx = m.pattern().as_usize();
                let p = m.start();
                let off1 = core_off[idx];
                if p < off1 {
                    continue;
                }
                let cand = p - off1;
                let masks_i = &masks[idx];
                let end = cand + masks_i.len();
                if end > seq.len() {
                    continue;
                }
                // 第二字面量段快速过滤（绝大多数候选在这里被拒绝）
                if let Some((off2, c2)) = &core2[idx] {
                    let s2 = cand + off2;
                    if seq[s2..s2 + c2.len()] != c2[..] {
                        continue;
                    }
                }
                // 逐位置掩码验证，等价于原正则的定长窗口匹配
                let window = &seq[cand..end];
                if window
                    .iter()
                    .zip(masks_i.iter())
                    .all(|(&b, &m)| m & BASE_BIT[b as usize] != 0)
                {
                    hits.push((idx, cand, enzyme.pattern_tag_lengths[idx]));
                }
            }
        }
    }
    // 归组排序：组内保持起点升序。
    hits.sort_by_key(|&(idx, s, _)| (idx, s));
    out.extend(hits.iter().map(|&(_, s, l)| (s, l)));
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
    /// 每个 tag 在其序列（contig）上的 bp 偏移（首个出现位置），与 `tags` 一一对应。
    /// 用于 anchor breadth/uniformity 置信度指标，也是向 TGT（Tag–Gap–Tag）
    /// 有序锚点模型对齐的第一步。v1 格式数据库无此字段，读取时为 None。
    #[serde(default)]
    pub tag_positions: Option<Vec<u32>>,
    /// 该 entry 对应序列（contig）的长度（bp）；0 表示未知（v1 格式数据库）。
    #[serde(default)]
    pub seq_len: u32,
}

// ---- .syldb 格式版本兼容 ----
//
// v2 起文件头带 magic，其后是 bincode 编码的 Vec<SyldbEntry>。
// v1（无 magic）是 bincode 编码的旧版 Vec<SyldbEntryV1>（无 tag_positions/seq_len），
// 读取时回退到旧 schema 并把新字段置为 None/0，保证既有数据库可直接加载
// （此时 profile 的 anchor 置信度列输出 NA）。
pub const SYLDB_MAGIC: &[u8; 8] = b"M2BDB\x00\x00\x02";

/// v1 格式的条目 schema（不含 tag_positions/seq_len），仅用于旧库回退读取。
/// （Serialize 仅供测试构造 v1 字节流。）
#[derive(Serialize, Deserialize)]
struct SyldbEntryV1 {
    sequence_id: String,
    tags: Vec<Hash>,
    tag_lengths: Vec<u8>,
    genome_source: String,
    tag_uniqueness: Option<Vec<bool>>,
    #[serde(default)]
    tag_seqs: Option<Vec<TagHash>>,
    #[serde(default)]
    enzyme: String,
}

impl From<SyldbEntryV1> for SyldbEntry {
    fn from(v1: SyldbEntryV1) -> Self {
        SyldbEntry {
            sequence_id: v1.sequence_id,
            tags: v1.tags,
            tag_lengths: v1.tag_lengths,
            genome_source: v1.genome_source,
            tag_uniqueness: v1.tag_uniqueness,
            tag_seqs: v1.tag_seqs,
            enzyme: v1.enzyme,
            tag_positions: None,
            seq_len: 0,
        }
    }
}

/// 写出 .syldb（v2：magic + bincode）。
pub fn write_syldb<W: Write>(mut writer: W, entries: &[SyldbEntry]) -> Result<()> {
    writer
        .write_all(SYLDB_MAGIC)
        .context("Failed to write syldb magic")?;
    bincode::serialize_into(writer, entries).context("Failed to serialize syldb data")
}

/// 读取 .syldb，自动识别 v2（带 magic）与 v1（旧版无 magic）格式。
pub fn read_syldb<R: Read>(mut reader: R) -> Result<Vec<SyldbEntry>> {
    let mut magic = [0u8; 8];
    if let Err(e) = reader.read_exact(&mut magic) {
        // 不足 8 字节的文件既不是 v2 也不可能是有效的 v1（空 Vec 也要 8 字节长度前缀）
        return Err(e).context("Failed to read syldb header (file too short)");
    }
    if &magic == SYLDB_MAGIC {
        use bincode::Options;
        let entries: Vec<SyldbEntry> = bincode::DefaultOptions::new()
            .with_fixint_encoding()
            .reject_trailing_bytes()
            .deserialize_from(reader)
            .context("Failed to deserialize v2 syldb payload")?;
        return Ok(entries);
    }
    // v1：把 magic 字节拼回流再按旧 schema 解析。reject_trailing_bytes 要求恰好
    // 消费整个流，误解析（非 syldb 的二进制文件）会迅速报错而不是错读。
    use bincode::Options;
    let chained = std::io::Cursor::new(magic.to_vec()).chain(reader);
    let v1: Vec<SyldbEntryV1> = bincode::DefaultOptions::new()
        .with_fixint_encoding()
        .reject_trailing_bytes()
        .deserialize_from(chained)
        .context("Failed to deserialize database file (unrecognized format; neither v2 magic nor legacy v1)")?;
    Ok(v1.into_iter().map(SyldbEntry::from).collect())
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
    // 每条记录复用的缓冲区，避免 per-record 分配
    let mut bufs = TagBufs::default();
    let mut tags: Vec<(Hash, TagHash, u8)> = Vec::with_capacity(64);

    while let Some(record) = reader.next() {
        if record.is_ok() {
            let record = record.expect(&format!("Invalid record for file {} ", input.display()));
            let seq = record.seq();
            let seq_id = String::from_utf8_lossy(record.id());

            stats.total_sequences += 1;
            stats.total_sequence_length += seq.len();

            // canonical tag 字节序列（用于写出 FASTA/FASTQ）
            extract_canonical_tags_into(&seq, enzyme, &mut bufs, &mut tags);

            // 按照sylph的去重模式（现在使用canonical tags）
            for (_, tag, _len) in tags.drain(..) {
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
    // 每条记录复用的缓冲区，避免 per-record 分配
    let mut bufs = TagBufs::default();
    let mut tags: Vec<(Hash, TagHash, u8)> = Vec::with_capacity(64);

    while let Some(record) = reader.next() {
        if record.is_ok() {
            let record = record.expect(&format!("Invalid record for file {} ", input.display()));
            let seq = record.seq();
            let seq_id = String::from_utf8_lossy(record.id());

            stats.total_sequences += 1;
            stats.total_sequence_length += seq.len();

            // canonical tag 字节序列（用于写出 FASTA/FASTQ）
            extract_canonical_tags_into(&seq, enzyme, &mut bufs, &mut tags);

            // 按照sylph的去重模式（现在使用canonical tags）
            for (_, tag, _len) in tags.drain(..) {
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

            // needletail 自动识别 gzip 压缩与 fasta/fastq 格式。
            // <2 字节的（空）文件旧版会得到 0 条记录，这里显式保持该行为。
            let mut stats = ExtractionStats::new();
            if get_file_size_optimized(&input_path)? < 2 {
                log_stats(stats, &enzyme);
                continue;
            }
            let mut reader = parse_fastx_file(&input_path)
                .context(format!("Failed to open reads file: {}", input_path.display()))?;

            // 单文件内按 chunk 并行提取：串行解析（gzip 解压本身无法并行），
            // 把一个 chunk 的 (id, seq) 拷出后由 rayon 并行做酶切扫描，
            // indexed par_iter + collect 保证结果顺序与输入一致，
            // 因此 fasta 与 .sylsp 输出与串行基线逐字节一致。
            const READ_CHUNK: usize = 2048;
            let mut chunk: Vec<(Vec<u8>, Vec<u8>)> = Vec::with_capacity(READ_CHUNK);
            'chunk_loop: loop {
                chunk.clear();
                while chunk.len() < READ_CHUNK {
                    match reader.next() {
                        Some(rec) => {
                            let rec = rec.context("Failed to read FASTQ record")?;
                            let seq = rec.seq();
                            stats.total_sequences += 1;
                            stats.total_sequence_length += seq.len();
                            chunk.push((fastq_id(rec.id()).to_vec(), seq.into_owned()));
                        }
                        None => {
                            if chunk.is_empty() {
                                break 'chunk_loop;
                            }
                            break;
                        }
                    }
                }

                // 每条 read 产出：fasta 片段字节 + 对应的 sylsp 条目
                let chunk_results: Vec<(Vec<u8>, Vec<SylspEntry>)> = chunk
                    .par_iter()
                    .map_init(TagBufs::default, |bufs, (id, seq)| {
                        let mut tags: Vec<(Hash, TagHash, u8)> = Vec::with_capacity(8);
                        extract_canonical_tags_into(seq, &enzyme, bufs, &mut tags);
                        let id_lossy = String::from_utf8_lossy(id);
                        let mut frag = Vec::new();
                        let mut entries = Vec::with_capacity(tags.len());
                        for (i, (h, tag, tag_len)) in tags.into_iter().enumerate() {
                            let entry_id = format!("{}_tag{}", id_lossy, i + 1);
                            let _ = writeln!(frag, ">{}\n{}", entry_id, String::from_utf8_lossy(&tag));
                            entries.push(SylspEntry {
                                sequence_id: entry_id,
                                tag: h,
                                tag_length: tag_len,
                                sample_source: file_stem.clone(),
                                enzyme: enzyme.name.clone(),
                            });
                        }
                        (frag, entries)
                    })
                    .collect();

                for (frag, entries) in chunk_results {
                    stats.total_tags += entries.len();
                    fa_writer.write_all(&frag).context("Failed to write FASTA record")?;
                    all_sylsp_entries.extend(entries);
                }
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
            
            write_syldb(combined_syldb_writer, &all_syldb_entries)
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
            
            write_syldb(combined_syldb_writer, &all_syldb_entries)
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

                let mut sylsp_entries = Vec::new();
                let mut stats = ExtractionStats::new();

                // <2 字节的（空）文件旧版会得到 0 条记录，保持该行为。
                if get_file_size_optimized(&input_path)? < 2 {
                    let mut global_stats = sample_stats.lock().unwrap();
                    global_stats.insert(file_stem.clone(), stats.clone());
                    log_stats(stats, &enzyme);
                    return Ok((file_stem, sylsp_entries));
                }
                // needletail 自动识别 gzip 与 fasta/fastq 格式
                let mut reader = parse_fastx_file(&input_path)
                    .context(format!("Failed to open sample file: {}", input_path.display()))?;

                // 每条 read 复用的缓冲区，避免 per-read 的 Vec/FxHashSet 分配
                let mut bufs = TagBufs::default();
                let mut tags: Vec<(Hash, u8)> = Vec::with_capacity(64);
                while let Some(rec) = reader.next() {
                    let rec = rec.context("Failed to read FASTQ record")?;
                    let seq = rec.seq();
                    stats.total_sequences += 1;
                    stats.total_sequence_length += seq.len();

                    extract_tag_hashes_into(&seq, &enzyme, &mut bufs, &mut tags);
                    stats.total_tags += tags.len();

                    let id = String::from_utf8_lossy(fastq_id(rec.id()));
                    for (i, (tag, tag_len)) in tags.iter().enumerate() {
                        sylsp_entries.push(SylspEntry {
                            sequence_id: format!("{}_tag{}", id, i + 1),
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

    // needletail 解析 FASTA（自动识别 gzip）；<2 字节的（空）文件按 0 条记录处理，
    // 与旧版 rust-bio 行为一致。seq() 已去除换行，等价于 rust-bio 的多行拼接。
    if get_file_size_optimized(input)? >= 2 {
        let mut reader = parse_fastx_file(input)
            .context(format!("Failed to open FASTA file: {}", input.display()))?;
        // 每条 contig 复用的缓冲区，避免 per-record 的 Vec/FxHashSet 分配
        let mut bufs = TagBufs::default();
        let mut tag_items: Vec<(Hash, TagHash, u8, u32)> = Vec::with_capacity(64);
        while let Some(rec) = reader.next() {
            let rec = rec.context("Failed to read FASTA record")?;
            let seq = rec.seq();
            stats.total_sequences += 1;
            stats.total_sequence_length += seq.len();

            // 提取 canonical tag 字节序列及其哈希；保留序列以支持 error-tolerant matching，
            // 保留 bp 位置以支持 anchor breadth/uniformity 置信度指标（TGT 对齐）。
            // 哈希在去重时已算好，直接复用，不再二次 hash_bytes。
            extract_canonical_tags_pos_into(&seq, enzyme, &mut bufs, &mut tag_items);
            stats.total_tags += tag_items.len();
            let mut tags = Vec::with_capacity(tag_items.len());
            let mut tag_lengths = Vec::with_capacity(tag_items.len());
            let mut tag_positions = Vec::with_capacity(tag_items.len());
            let mut tag_seqs = if store_tag_seqs {
                Some(Vec::with_capacity(tag_items.len()))
            } else {
                None
            };
            for (h, tag, tag_len, pos) in tag_items.drain(..) {
                tags.push(h);
                tag_lengths.push(tag_len);
                tag_positions.push(pos);
                if let Some(seqs) = tag_seqs.as_mut() {
                    seqs.push(tag);
                }
            }

            // 创建 syldb 条目
            let entry = SyldbEntry {
                sequence_id: String::from_utf8_lossy(fasta_id(rec.id())).into_owned(),
                tags,
                tag_lengths,
                genome_source: input.to_string_lossy().to_string(),
                tag_uniqueness: None, // 初始时未标记，将由mark命令处理
                tag_seqs,
                enzyme: enzyme_name.clone(),
                tag_positions: Some(tag_positions),
                seq_len: seq.len() as u32,
            };
            syldb_entries.push(entry);
        }
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
    // needletail 自动识别 gzip；<2 字节的（空）文件按 0 条记录处理，与旧版一致。
    let mut stats = ExtractionStats::new();
    let mut entries = Vec::new();

    if get_file_size_optimized(Path::new(input1))? < 2 || get_file_size_optimized(Path::new(input2))? < 2 {
        log_stats(stats, enzyme);
        return Ok(entries);
    }
    let mut reader1 = parse_fastx_file(Path::new(input1))
        .context(format!("Failed to open first pair file: {}", input1))?;
    let mut reader2 = parse_fastx_file(Path::new(input2))
        .context(format!("Failed to open second pair file: {}", input2))?;

    // 每条 read 复用的缓冲区，避免 per-read 的 Vec/FxHashSet 分配
    let mut bufs = TagBufs::default();
    let mut tags1: Vec<(Hash, u8)> = Vec::with_capacity(8);
    let mut tags2: Vec<(Hash, u8)> = Vec::with_capacity(8);

    loop {
        let record1 = match reader1.next() {
            Some(Ok(r)) => r,
            Some(Err(e)) => return Err(anyhow::anyhow!("Error reading first pair: {}", e)),
            None => break,
        };

        let record2 = match reader2.next() {
            Some(Ok(r)) => r,
            Some(Err(e)) => return Err(anyhow::anyhow!("Error reading second pair: {}", e)),
            None => break,
        };

        let seq1 = record1.seq();
        let seq2 = record2.seq();
        stats.total_sequences += 1;
        stats.total_sequence_length += seq1.len() + seq2.len();

        // 每条 read 内 tag 已按哈希去重；read id + tag index 天然唯一，
        // 原先的 seen_pairs 集合永远不会命中重复，纯属浪费，已移除。
        extract_tag_hashes_into(&seq1, enzyme, &mut bufs, &mut tags1);
        extract_tag_hashes_into(&seq2, enzyme, &mut bufs, &mut tags2);
        stats.total_tags += tags1.len() + tags2.len();

        let id1 = String::from_utf8_lossy(fastq_id(record1.id())).into_owned();
        for (i, (tag, tag_len)) in tags1.iter().enumerate() {
            entries.push((format!("{}_{}", id1, i + 1), *tag, *tag_len, sample_source.to_string()));
        }
        let id2 = String::from_utf8_lossy(fastq_id(record2.id())).into_owned();
        for (i, (tag, tag_len)) in tags2.iter().enumerate() {
            entries.push((format!("{}_{}", id2, i + 1), *tag, *tag_len, sample_source.to_string()));
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

    /// 旧版逐 pattern 扫描实现，仅作为单趟扫描等价性测试的对照。
    fn find_all_tag_positions_per_pattern(seq: &[u8], enzyme: &EnzymeSpec) -> Vec<(usize, usize)> {
        let mut out = Vec::new();
        for (pattern, &tag_len) in enzyme.patterns.iter().zip(&enzyme.pattern_tag_lengths) {
            let mut start = 0usize;
            while start <= seq.len() {
                match pattern.find_at(seq, start) {
                    Some(m) => {
                        let mstart = m.start();
                        out.push((mstart, tag_len));
                        start = mstart + 1;
                    }
                    None => break,
                }
            }
        }
        out
    }

    #[test]
    fn test_single_pass_scan_matches_per_pattern() {
        // 确定性伪随机序列（xorshift64*），覆盖多种长度
        let mut state = 0x9E3779B97F4A7C15u64;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };

        // 说明：密集位点种子在下方通过对长随机序列的真实命中收集，
        // 用于人为制造重叠/相邻命中，压力测试同起点与交错命中的归组顺序。

        for enzyme_name in ["BcgI", "all"] {
            let enzyme = EnzymeSpec::new(enzyme_name).unwrap();

            // 1) 纯随机序列
            for len in [0usize, 10, 33, 100, 1000, 5000] {
                for _ in 0..20 {
                    let seq: Vec<u8> = (0..len).map(|_| b"ACGT"[(next() % 4) as usize]).collect();
                    assert_eq!(
                        find_all_tag_positions_per_pattern(&seq, &enzyme),
                        find_all_tag_positions(&seq, &enzyme),
                        "enzyme={} len={}",
                        enzyme_name,
                        len
                    );
                }
            }

            // 2) 先扫一段长随机序列收集真实命中，再把命中片段密集拼接到新序列里，
            //    人为制造大量（重叠/相邻）位点
            let big: Vec<u8> = (0..20000).map(|_| b"ACGT"[(next() % 4) as usize]).collect();
            let seeds = find_all_tag_positions_per_pattern(&big, &enzyme);
            if !seeds.is_empty() {
                let mut dense = Vec::new();
                for &(s, l) in seeds.iter().take(200) {
                    dense.extend_from_slice(&big[s..s + l]);
                    // 随机回退 0..l-1 bp，使下一段与上一段（部分）重叠
                    let back = (next() as usize) % l;
                    dense.truncate(dense.len() - back.min(dense.len()));
                }
                assert_eq!(
                    find_all_tag_positions_per_pattern(&dense, &enzyme),
                    find_all_tag_positions(&dense, &enzyme),
                    "enzyme={} dense",
                    enzyme_name
                );
            }
        }
    }

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

    // ---- 多成员 gzip 回归测试 ----
    // 历史 bug：旧实现用 flate2::GzDecoder 只读 gzip 第一个成员，
    // `cat R1.gz R2.gz` 拼接的输入会静默丢掉后半部分 reads。
    // 换成 needletail 后已修复，以下测试防止回归。

    /// 将 data 单独压缩成一个 gzip 成员的字节流
    fn gzip_member(data: &[u8]) -> Vec<u8> {
        use flate2::write::GzEncoder;
        let mut enc = GzEncoder::new(Vec::new(), flate2::Compression::fast());
        std::io::Write::write_all(&mut enc, data).unwrap();
        enc.finish().unwrap()
    }

    /// 第 i 条含单个 BcgI 位点的 read：10bp + CGA + 6bp + TGC + 10bp = 32bp，
    /// 可变区用 i 的四进制编码填充，保证每条 read 的 tag 唯一（避免全局去重干扰计数）。
    fn bcgi_read(i: usize) -> String {
        let b4 = |mut v: usize, len: usize| -> String {
            let mut s = String::with_capacity(len);
            for _ in 0..len {
                s.push(b"ACGT"[v % 4] as char);
                v /= 4;
            }
            s
        };
        format!("{}CGA{}TGC{}", b4(i, 10), b4(i, 6), b4(i.wrapping_add(12345), 10))
    }

    fn fastq_records(n: usize, offset: usize) -> String {
        let mut s = String::new();
        for i in 0..n {
            let read = bcgi_read(offset + i);
            s.push_str(&format!("@r{}\n{}\n+\n{}\n", offset + i, read, "~".repeat(read.len())));
        }
        s
    }

    /// 统计 extract 输出 FASTA 的 tag 条数（每条 tag 一行 header）
    fn count_fa_tags(path: &Path) -> usize {
        let content = std::fs::read_to_string(path).unwrap();
        content.lines().filter(|l| l.starts_with('>')).count()
    }

    #[test]
    fn test_multi_member_gzip_fastq_reads_all_members() {
        let dir = std::env::temp_dir();
        let n = 50;
        let member1 = fastq_records(n, 0);
        let member2 = fastq_records(n, 1000);
        let enzyme = EnzymeSpec::new("BcgI").unwrap();

        // 对照：两个单成员文件分别处理
        let m1_path = dir.join("m2b_test_mm_m1.fq.gz");
        let m2_path = dir.join("m2b_test_mm_m2.fq.gz");
        std::fs::write(&m1_path, gzip_member(member1.as_bytes())).unwrap();
        std::fs::write(&m2_path, gzip_member(member2.as_bytes())).unwrap();
        let out1 = dir.join("m2b_test_mm_out1.fa");
        let out2 = dir.join("m2b_test_mm_out2.fa");
        process_fastq(&m1_path, &out1, &enzyme, "fa", false).unwrap();
        process_fastq(&m2_path, &out2, &enzyme, "fa", false).unwrap();
        let count1 = count_fa_tags(&out1);
        let count2 = count_fa_tags(&out2);
        assert_eq!(count1, n, "member1 应提取 {} 条 tag", n);
        assert_eq!(count2, n, "member2 应提取 {} 条 tag", n);

        // 被测：cat m1.gz m2.gz 的多成员文件必须读到两个成员
        let mut cat_bytes = std::fs::read(&m1_path).unwrap();
        cat_bytes.extend_from_slice(&std::fs::read(&m2_path).unwrap());
        let cat_path = dir.join("m2b_test_mm_cat.fq.gz");
        std::fs::write(&cat_path, &cat_bytes).unwrap();
        let out_cat = dir.join("m2b_test_mm_outcat.fa");
        process_fastq(&cat_path, &out_cat, &enzyme, "fa", false).unwrap();
        let count_cat = count_fa_tags(&out_cat);
        assert_eq!(
            count_cat,
            count1 + count2,
            "多成员 gzip 只读到 {} 条 tag（期望 {}），疑似只读了第一个成员",
            count_cat,
            count1 + count2
        );

        for p in [&m1_path, &m2_path, &out1, &out2, &cat_path, &out_cat] {
            std::fs::remove_file(p).ok();
        }
    }

    #[test]
    fn test_multi_member_gzip_fasta_reads_all_members() {
        let dir = std::env::temp_dir();
        let enzyme = EnzymeSpec::new("BcgI").unwrap();
        // 每个成员一条“基因组”，各含 20 个 BcgI 位点
        let genome = |offset: usize| -> String {
            let mut s = String::from(">genome\n");
            for i in 0..20 {
                s.push_str(&bcgi_read(offset + i));
            }
            s.push('\n');
            s
        };
        let m1 = gzip_member(genome(0).as_bytes());
        let m2 = gzip_member(genome(100000).as_bytes());
        let mut cat_bytes = m1.clone();
        cat_bytes.extend_from_slice(&m2);
        let cat_path = dir.join("m2b_test_mm_genomes.fa.gz");
        std::fs::write(&cat_path, &cat_bytes).unwrap();
        let out_base = dir.join("m2b_test_mm_genomes");

        let sketches = process_fasta_to_syldb(&cat_path, &out_base, &enzyme, "fa", false, false).unwrap();
        assert_eq!(sketches.len(), 2, "多成员 gzip FASTA 应读出 2 条基因组");
        // 每条基因组至少 20 个 tag（read 拼接处可能偶然形成额外位点，故不取等号）
        assert!(
            sketches.iter().all(|g| g.tags.len() >= 20),
            "每条基因组应至少有 20 个 tag，实际 {:?}",
            sketches.iter().map(|g| g.tags.len()).collect::<Vec<_>>()
        );

        std::fs::remove_file(&cat_path).ok();
        std::fs::remove_file(out_base.with_extension("syldb")).ok();
    }

    // ---- syldb v2 格式（位置索引 + magic）与 v1 回退 ----

    #[test]
    fn test_syldb_write_read_roundtrip_v2() {
        let entry = SyldbEntry {
            sequence_id: "g1".to_string(),
            tags: vec![11, 22, 33],
            tag_lengths: vec![32, 32, 32],
            genome_source: "g1.fasta".to_string(),
            tag_uniqueness: None,
            tag_seqs: Some(vec![b"AAAA".to_vec(), b"CCCC".to_vec(), b"GGGG".to_vec()]),
            enzyme: "BcgI".to_string(),
            tag_positions: Some(vec![100, 2000, 5000]),
            seq_len: 8000,
        };
        let mut buf = Vec::new();
        write_syldb(&mut buf, std::slice::from_ref(&entry)).unwrap();
        assert!(buf.starts_with(SYLDB_MAGIC), "v2 文件应以 magic 开头");

        let back = read_syldb(std::io::Cursor::new(&buf)).unwrap();
        assert_eq!(back.len(), 1);
        assert_eq!(back[0].tags, entry.tags);
        assert_eq!(back[0].tag_positions, entry.tag_positions);
        assert_eq!(back[0].seq_len, entry.seq_len);
        assert_eq!(back[0].enzyme, entry.enzyme);
    }

    #[test]
    fn test_syldb_read_v1_fallback() {
        // v1 字节流（无 magic、无 tag_positions/seq_len）必须可读，新字段置 None/0
        let v1 = SyldbEntryV1 {
            sequence_id: "g1".to_string(),
            tags: vec![7, 8],
            tag_lengths: vec![32, 32],
            genome_source: "g1.fasta".to_string(),
            tag_uniqueness: None,
            tag_seqs: None,
            enzyme: "BcgI".to_string(),
        };
        let buf = bincode::serialize(&vec![v1]).unwrap();
        let back = read_syldb(std::io::Cursor::new(&buf)).unwrap();
        assert_eq!(back.len(), 1);
        assert_eq!(back[0].tags, vec![7, 8]);
        assert!(back[0].tag_positions.is_none());
        assert_eq!(back[0].seq_len, 0);
    }

    #[test]
    fn test_process_fasta_to_syldb_stores_positions() {
        let dir = std::env::temp_dir();
        let fa_path = dir.join("m2b_test_pos.fa");
        // 两条 contig，各含若干 BcgI 位点
        let mut contig1 = String::new();
        for i in 0..5 {
            contig1.push_str(&bcgi_read(i));
        }
        let mut contig2 = String::new();
        for i in 0..3 {
            contig2.push_str(&bcgi_read(100 + i));
        }
        std::fs::write(
            &fa_path,
            format!(">c1\n{}\n>c2\n{}\n", contig1, contig2),
        )
        .unwrap();
        let out_base = dir.join("m2b_test_pos");
        let enzyme = EnzymeSpec::new("BcgI").unwrap();
        let entries = process_fasta_to_syldb(&fa_path, &out_base, &enzyme, "fa", false, true).unwrap();
        assert_eq!(entries.len(), 2);
        for (e, expected_len) in entries.iter().zip([contig1.len(), contig2.len()]) {
            let pos = e.tag_positions.as_ref().expect("v2 应存储位置");
            assert_eq!(pos.len(), e.tags.len(), "位置与 tag 一一对应");
            assert_eq!(e.seq_len as usize, expected_len, "seq_len 应为 contig 长度");
            assert!(
                pos.iter().all(|&p| (p as usize) < expected_len),
                "位置必须落在 contig 内"
            );
        }
        std::fs::remove_file(&fa_path).ok();
    }
}
