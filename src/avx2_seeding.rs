use std::arch::x86_64::*;
use crate::constants::*;

#[inline]
#[target_feature(enable = "avx2")]
pub unsafe fn mm_hash256(kmer: __m256i) -> __m256i {
    let mut key = kmer;
    let s1 = _mm256_slli_epi64(key, 21);
    key = _mm256_add_epi64(key, s1);
    key = _mm256_xor_si256(key, _mm256_cmpeq_epi64(key, key));

    key = _mm256_xor_si256(key, _mm256_srli_epi64(key, 24));
    let s2 = _mm256_slli_epi64(key, 3);
    let s3 = _mm256_slli_epi64(key, 8);

    key = _mm256_add_epi64(key, s2);
    key = _mm256_add_epi64(key, s3);
    key = _mm256_xor_si256(key, _mm256_srli_epi64(key, 14));
    let s4 = _mm256_slli_epi64(key, 2);
    let s5 = _mm256_slli_epi64(key, 4);
    key = _mm256_add_epi64(key, s4);
    key = _mm256_add_epi64(key, s5);
    key = _mm256_xor_si256(key, _mm256_srli_epi64(key, 28));

    let s6 = _mm256_slli_epi64(key, 31);
    key = _mm256_add_epi64(key, s6);

    return key;
}

// ---- 与 sketch::extract_kmers / extract_kmers_positions 逐位等价的 AVX2 版本 ----
// 滚动 k-mer 更新保持标量逻辑（含无效碱基跳过语义），仅把 mm_hash64 与阈值比较
// 按 4 通道向量化。原 extract_markers_avx2（从 sylph 移植）语义与 sketch 热循环
// 不一致（仅支持 k=21/31、不跳过 N、丢弃尾部余量、分块导致输出顺序不同），已删除。

// 无符号 h < threshold 的 4 通道比较（AVX2 只有有符号 cmpgt，先翻符号位）
#[inline]
#[target_feature(enable = "avx2")]
unsafe fn hash_lt_threshold(hash: __m256i, threshold: u64) -> [i64; 4] {
    let sign = _mm256_set1_epi64x(i64::MIN);
    let thr = _mm256_xor_si256(_mm256_set1_epi64x(threshold as i64), sign);
    let m = _mm256_cmpgt_epi64(thr, _mm256_xor_si256(hash, sign));
    let mut mv = [0i64; 4];
    _mm256_storeu_si256(mv.as_mut_ptr() as *mut __m256i, m);
    mv
}

#[target_feature(enable = "avx2")]
pub unsafe fn extract_kmers_avx2(string: &[u8], kmer_vec: &mut Vec<Hash>, c: usize, k: usize) {
    if string.len() < k {
        return;
    }

    let mut rolling_kmer_f: u64 = 0;
    let mut rolling_kmer_r: u64 = 0;

    let reverse_shift_dist = 2 * (k - 1);
    let mask = u64::MAX >> (std::mem::size_of::<u64>() * 8 - 2 * k);
    let rev_mask = !(3 << (2 * k - 2));
    let len = string.len();
    let threshold = u64::MAX / (c as u64);

    // 初始化前 k-1 个核苷酸
    for i in 0..k - 1 {
        let nuc_f = crate::sketch::BYTE_TO_SEQ[string[i] as usize] as u64;
        if nuc_f >= 4 {
            return; // 跳过包含无效核苷酸的序列
        }
        let nuc_r = 3 - nuc_f;
        rolling_kmer_f <<= 2;
        rolling_kmer_f |= nuc_f;
        rolling_kmer_r >>= 2;
        rolling_kmer_r |= nuc_r << reverse_shift_dist;
    }

    let mut buf = [0u64; 4];
    let mut nbuf = 0usize;

    // 滑动窗口提取k-mers
    for i in k - 1..len {
        let nuc_byte = string[i] as usize;
        let nuc_f = crate::sketch::BYTE_TO_SEQ[nuc_byte] as u64;
        if nuc_f >= 4 {
            continue; // 跳过无效核苷酸
        }
        let nuc_r = 3 - nuc_f;

        rolling_kmer_f <<= 2;
        rolling_kmer_f |= nuc_f;
        rolling_kmer_f &= mask;

        rolling_kmer_r >>= 2;
        rolling_kmer_r &= rev_mask;
        rolling_kmer_r |= nuc_r << reverse_shift_dist;

        // 选择canonical k-mer
        let canonical_kmer = if rolling_kmer_f < rolling_kmer_r {
            rolling_kmer_f
        } else {
            rolling_kmer_r
        };

        buf[nbuf] = canonical_kmer;
        nbuf += 1;
        if nbuf == 4 {
            let h = mm_hash256(_mm256_loadu_si256(buf.as_ptr() as *const __m256i));
            let mv = hash_lt_threshold(h, threshold);
            let mut hv = [0u64; 4];
            _mm256_storeu_si256(hv.as_mut_ptr() as *mut __m256i, h);
            for lane in 0..4 {
                if mv[lane] != 0 {
                    kmer_vec.push(hv[lane]);
                }
            }
            nbuf = 0;
        }
    }

    // 尾部余量用标量处理
    for &km in &buf[..nbuf] {
        let hash_value = crate::sketch::mm_hash64(km);
        if hash_value < threshold {
            kmer_vec.push(hash_value);
        }
    }
}

#[target_feature(enable = "avx2")]
pub unsafe fn extract_kmers_positions_avx2(
    string: &[u8],
    kmer_vec: &mut Vec<(usize, usize, u64)>,
    c: usize,
    k: usize,
    contig_number: usize,
) {
    if string.len() < k {
        return;
    }

    let mut rolling_kmer_f: u64 = 0;
    let mut rolling_kmer_r: u64 = 0;

    let reverse_shift_dist = 2 * (k - 1);
    let mask = u64::MAX >> (std::mem::size_of::<u64>() * 8 - 2 * k);
    let rev_mask = !(3 << (2 * k - 2));
    let len = string.len();
    let threshold = u64::MAX / (c as u64);

    // 初始化前 k-1 个核苷酸
    for i in 0..k - 1 {
        let nuc_f = crate::sketch::BYTE_TO_SEQ[string[i] as usize] as u64;
        if nuc_f >= 4 {
            return; // 跳过包含无效核苷酸的序列
        }
        let nuc_r = 3 - nuc_f;
        rolling_kmer_f <<= 2;
        rolling_kmer_f |= nuc_f;
        rolling_kmer_r >>= 2;
        rolling_kmer_r |= nuc_r << reverse_shift_dist;
    }

    let mut buf = [0u64; 4];
    let mut pos_buf = [0usize; 4];
    let mut nbuf = 0usize;

    // 滑动窗口提取k-mers
    for i in k - 1..len {
        let nuc_byte = string[i] as usize;
        let nuc_f = crate::sketch::BYTE_TO_SEQ[nuc_byte] as u64;
        if nuc_f >= 4 {
            continue; // 跳过无效核苷酸
        }
        let nuc_r = 3 - nuc_f;

        rolling_kmer_f <<= 2;
        rolling_kmer_f |= nuc_f;
        rolling_kmer_f &= mask;

        rolling_kmer_r >>= 2;
        rolling_kmer_r &= rev_mask;
        rolling_kmer_r |= nuc_r << reverse_shift_dist;

        // 选择canonical k-mer
        let canonical_kmer = if rolling_kmer_f < rolling_kmer_r {
            rolling_kmer_f
        } else {
            rolling_kmer_r
        };

        buf[nbuf] = canonical_kmer;
        pos_buf[nbuf] = i;
        nbuf += 1;
        if nbuf == 4 {
            let h = mm_hash256(_mm256_loadu_si256(buf.as_ptr() as *const __m256i));
            let mv = hash_lt_threshold(h, threshold);
            let mut hv = [0u64; 4];
            _mm256_storeu_si256(hv.as_mut_ptr() as *mut __m256i, h);
            for lane in 0..4 {
                if mv[lane] != 0 {
                    kmer_vec.push((contig_number, pos_buf[lane], hv[lane]));
                }
            }
            nbuf = 0;
        }
    }

    // 尾部余量用标量处理
    for lane in 0..nbuf {
        let hash_value = crate::sketch::mm_hash64(buf[lane]);
        if hash_value < threshold {
            kmer_vec.push((contig_number, pos_buf[lane], hash_value));
        }
    }
}
