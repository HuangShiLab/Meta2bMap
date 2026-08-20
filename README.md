# meta2bmap

**Ultrafast coverage-adjusted ANI queries and taxonomic profiling for 2bRAD and shotgun metagenomic data.**

meta2bmap is a Rust command-line tool that profiles metagenomic samples by matching
restriction-enzyme (2bRAD) tags extracted from sequencing reads against tags extracted
from reference genomes. It reports per-genome containment ANI with a coverage-adjustment
model (a negative-binomial per-locus model that corrects the naive containment ANI at low
coverage), together with taxonomic and sequence abundances.

It works on two kinds of data:

- **2bRAD / type-IIB restriction-enzyme data** — reads are (or contain) enzyme tags;
  meta2bmap extracts the tags directly.
- **Shotgun data** — the same enzyme-tag scheme is applied in silico to shotgun reads,
  yielding a sparse but fast marker set.

meta2bmap is part of the 2bRAD software family (with Fast2bRAD-M and Strain2bScan) and is
a derivative of [sylph](https://github.com/bluenote-1577/sylph) (see
[Relationship to sylph](#relationship-to-sylph) below).

Repository: <https://github.com/HuangShiLab/Meta2bMap>

## Installation

meta2bmap is built from source with the Rust toolchain (Rust 1.70+).

### Supported platforms

- macOS (10.15+)
- Linux: Ubuntu 16.04+ (incl. WSL), CentOS 7+, other distributions with glibc 2.17+

### Install Rust

```bash
# macOS: install Xcode command line tools first
xcode-select --install

# All platforms: install Rust via rustup
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source ~/.cargo/env

# Ubuntu/Debian (incl. WSL) also need:
sudo apt update && sudo apt install build-essential curl

# CentOS/RHEL equivalents:
sudo yum groupinstall "Development Tools" && sudo yum install curl
```

### Build

```bash
git clone https://github.com/HuangShiLab/Meta2bMap.git
cd Meta2bMap
cargo build --release
```

The binary is produced at `target/release/meta2bmap`.

## Quickstart

```bash
# 1. Extract 2bRAD tags from reference genomes AND a sample (one command)
meta2bmap extract -t 8 \
  -g genome1.fa genome2.fa \
  -r sample.fq \
  -o out -d out -n myrun -e BcgI

# -> out/myrun.syldb  (genome database)
# -> out/myrun.sylsp  (sample tags)

# 2. Profile: species-level abundances + coverage-adjusted ANI
meta2bmap profile \
  --db-file out/myrun.syldb \
  --sample-file out/myrun.sylsp \
  --minimum-ani 90 \
  --threads 8 \
  --tsv-name out/myrun_abundance.tsv > out/myrun_profile.txt
```

The genome composition table goes to **stdout** (or `--out-file-name`); the abundance
matrix goes to `--tsv-name`. Note that `profile` takes `--threads` (long form only),
while `extract` uses `-t/--threads`.

Other useful `extract` input modes:

```bash
# paired-end reads
meta2bmap extract -t 8 -1 a_1.fq b_1.fq -2 a_2.fq b_2.fq -d out -n paired

# batch mode from list files
meta2bmap extract -t 8 -k genome_list.txt -s reads_list.txt -o out -d out -n batch
```

`-k/--genome-list` accepts a two-column TSV (`genome_id<TAB>fasta_path`) or one path per
line; `-s/--sample-list` is one fastq path per line. Use `--l1`/`--l2` for paired list
files.

Useful `extract` options:

- `--no-tag-seqs`: do not store tag sequences in the `.syldb` (saves ~4-5x memory/disk
  for large databases; `--mismatch 1` will be unavailable for that database).
- `--max-ram <GB>`: cap RAM usage (default: 75% of total system memory; the extractor
  throttles itself to stay under this limit).

## Large databases: chunked build + merge

For very large reference collections, build per-chunk databases in parallel and merge
them (all chunks must use the same enzyme set):

```bash
meta2bmap extract -t 16 -k chunk1_list.txt -o chunk1 -d chunk1 -n chunk1 -e CjeI,BslFI
meta2bmap extract -t 16 -k chunk2_list.txt -o chunk2 -d chunk2 -n chunk2 -e CjeI,BslFI
meta2bmap merge chunk1/chunk1.syldb chunk2/chunk2.syldb -o full.syldb
```

`merge` refuses to combine databases built with different enzyme sets and warns if tag
sequences are stored inconsistently across inputs.

## Subcommands

| Subcommand | Purpose |
|------------|---------|
| `extract`  | Extract enzyme tags from FASTA genomes / FASTQ reads into `.syldb` / `.sylsp` |
| `merge`    | Merge multiple `.syldb` chunk files into one database (for parallel batch builds) |
| `profile`  | Species-level taxonomic profiling with abundances and coverage-adjusted ANI |
| `query`    | Coverage-adjusted containment-ANI queries between databases and samples |
| `sketch`   | sylph-style k-mer sampling sketching (alternative to tag extraction) |
| `inspect`  | Inspect extracted `.syldb`/`.sylsp` files (YAML summary / tag count matrix) |
| `view`     | View sketched `.db`/`.sp` files produced by `sketch` |
| `mark`     | Mark unique (taxa-specific) tags in a `.syldb` file |

Run `meta2bmap <subcommand> -h` for the full option list.

### Key `profile` options

- `--db-file`, `--sample-file` (required): pre-extracted `.syldb` / `.sylsp` files.
- `--minimum-ani <0-100>`: minimum adjusted ANI to report (default 95; values below 95
  give less accurate results).
- `--mismatch 0|1`: allow up to 1 Hamming mismatch when matching sample tags to database
  tags (more robust to sequencing error / strain variation; default 0).
- `--threads <N>`: threads (default 1; **no short form**).
- `--tsv-name <FILE>`: abundance matrix output (default `abundance_matrix.tsv`).
- `--taxonomy-file <FILE>`: GTDB-style taxonomy file for species-level aggregation.
- `--gscore-threshold <F>`: minimum G-score `sqrt(reads_count * tag_count)` for species
  filtering (default 10).
- `--min-eff-coverage <F>`: minimum effective coverage to retain a genome (default 0 =
  no filter).
- `--min-shared-tags <N>`: minimum shared tags required to report a genome (default 100;
  lower values increase sensitivity on large dense databases).
- `--min-tags-genome <N>`: exclude genomes with fewer than N matched tags up front
  (default 50).
- `--sketch-mode`: treat `--db-file`/`--sample-file` as k-mer sketches (`.db`/`.sp`)
  produced by `sketch` instead of 2bRAD tag files, and run them through the same
  containment/ANI engine (requires `--mismatch 0`).
- `--read-error-rate <F>` / `--no-error-correction`: in `--mismatch 1` mode, override or
  disable the error-aware ANI inversion (by default the per-base read error rate is
  estimated from the top candidate genomes).
- `--allow-enzyme-mismatch`: proceed even if the sample was extracted with a different
  enzyme set than the database (default: hard error).
- `--reassign-protection`: opt-in guard against winner-take-all over-stripping (default
  OFF; on large dense databases it rescues strain-cluster losers alongside true
  positives, so it costs precision).

`query` takes pre-extracted files as positional arguments
(`meta2bmap query *.syldb *.sylsp > results.tsv`) and shares most algorithm options;
its output header includes `ANI(%)`, `Eff_cov`, `ANI_5-95%`, `Median_cov`, `Mean_cov`,
`Containment`, and `Naive_ANI` per genome–sample pair.

## Interpreting the profile output

The genome composition table (stdout) has one row per detected genome:

| Column | Meaning |
|--------|---------|
| `Genome_ID` / `Sample_ID` | Reference genome and sample identifiers |
| `ANI(%)` | Coverage-adjusted ANI estimate |
| `Tax_Abund(%)` | Taxonomic abundance (relative abundance of the taxon) |
| `Seq_Abund(%)` | Sequence abundance (fraction of reads/tags assigned) |
| `Common_Tags` / `Total_Tags` | Tags shared with the sample / total tags in the genome |
| `Eff_cov` | Effective per-locus coverage estimated by the coverage model |

The abundance matrix written to `--tsv-name` contains one row per genome × sample with
the same abundance estimates in TSV form.

## k-mer sketch mode (shotgun)

As an alternative to enzyme tags, `sketch` builds sylph-style sampled k-mer sketches
(`.db` / `.sp`), which `profile --sketch-mode` and `query --sketch-mode` run through the
same containment / coverage-adjusted ANI engine:

```bash
meta2bmap sketch -t 8 -g refs/*.fa -o out --db-out-name refdb
meta2bmap sketch -t 8 -r sample.fq -d out
meta2bmap profile --sketch-mode --db-file out/refdb.db \
  --sample-file out/sample.sp --minimum-ani 90 --threads 8
```

Key `sketch` options: `-c/--c-value` subsampling rate (default 200), `--k-size`
(default 31), `--min-spacing` between selected k-mers (default 30). For benchmarking
against sylph, `sketch --sylph-faithful` builds sylph-equivalent sketches (deduplicate
repeated genomic k-mers instead of deleting them, and force `--min-spacing 0`).

## Multi-enzyme usage

`-e/--enzyme` accepts a single enzyme name, a comma-separated list, or `all`:

```bash
meta2bmap extract -t 8 -g refs/*.fa -r sample.fq -o out -d out -n multi -e BcgI,BslFI,CspCI,AloI
meta2bmap extract -t 8 -g refs/*.fa -r sample.fq -o out -d out -n all16  -e all
```

16 type-IIB enzymes are built in, with tag lengths stored per tag in the database:
AlfI (32 bp), AloI (27), BaeI (28), BcgI (32), BplI (27), BsaXI (27), BslFI (25),
Bsp24I (27), CjeI (28), CjePI (27), CspCI (33), FalI (27), HaeIV (27), Hin4I (27),
PpiI (27), PsrI (27). Using several enzymes increases the number of markers per genome
and improves accuracy at the cost of larger databases and longer extraction time.

## File formats

- `.syldb` — serialized genome database: per-genome tag hashes with **per-tag lengths**
  (enabling joint multi-enzyme ANI estimation) and the enzyme set recorded in each
  entry, so the file is self-describing.
- `.sylsp` — serialized sample tags with per-tag lengths.
- `.db` / `.sp` — legacy / `sketch`-mode formats (inspected with `view`).

**The format has changed over development: re-extract any `.syldb`/`.sylsp` files
produced by older builds.**

## Relationship to sylph

meta2bmap is a derivative work of [sylph](https://github.com/bluenote-1577/sylph) by
Jim Shaw, reusing its containment-ANI framework and parts of its code (e.g. the AVX2
k-mer hashing), with modifications by the HuangShi Lab for restriction-enzyme tag
extraction, multi-enzyme joint ANI, and the coverage-adjustment model. sylph is
distributed under the MIT license; meta2bmap preserves sylph's copyright notice and
license terms in [LICENSE](LICENSE). If you use meta2bmap, please also credit sylph.

## Citation

A manuscript describing meta2bmap is in preparation. Until then, please cite this
repository and the sylph paper (Shaw & Yu, *Nature Biotechnology*, 2024).
