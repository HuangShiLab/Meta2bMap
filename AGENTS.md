# Meta2bMap — Agent Guide

This guide is intended for AI coding agents that need to work on the Meta2bMap codebase. It assumes no prior knowledge of the project.

## Project Overview

Meta2bMap is a Rust command-line tool for **ultrafast genome ANI queries and taxonomic profiling of metagenomic shotgun samples**. It is built around two main ideas:

1. **2bRAD tag extraction**: restriction-enzyme-based extraction of short tags from reference genomes and sequencing reads.
2. **k-mer sketching**: a sylph-style sampled k-mer mode for ANI estimation.

The tool currently supports three categories of analysis (as described in `README.md`):

- taxonomic profiling;
- strain and functional profiling;
- statistical analysis of containment/ANI results.

The repository is located at `/Users/macstudio/meta2bmap`. It is a single Cargo workspace with one binary crate (`meta2bmap`).

## Technology Stack

- **Language**: Rust (edition 2021).
- **Minimum Rust version**: 1.70.0 (stated in `README.md`).
- **Build system**: Cargo.
- **CLI parsing**: `clap` with derive macros.
- **Parallelism**: `rayon` for data parallelism.
- **Allocator**: `tikv-jemallocator` is set as the global allocator in `src/main.rs`.
- **Sequence I/O**: `bio::io::{fasta,fastq}` and `needletail::parse_fastx_file`.
- **Compression**: `flate2` for gzip handling.
- **Serialization**: `bincode` and `serde` for `.db`/`.sp`/`.syldb`/`.sylsp` files.
- **Logging**: `log`, `env_logger`, and `simple_logger`.
- **Hash maps/sets**: `fxhash::{FxHashMap,FxHashSet}`.
- **Approximate set membership**: `scalable_cuckoo_filter` (used in sketching for deduplication).
- **System/memory info**: `sysinfo`, `memory-stats`.
- **Other utilities**: `anyhow`, `regex` (including `regex::bytes`), `itertools`, `glob`.

## Build Commands

Build the release binary:

```bash
cd /Users/macstudio/meta2bmap
cargo build --release
```

The binary is produced at:

```
target/release/meta2bmap
```

Build for development (includes debug symbols; release profile also sets `debug = true`):

```bash
cargo build
```

Run tests:

```bash
cargo test
```

Current test coverage is modest: unit tests live in `#[cfg(test)] mod tests` blocks inside `src/contain.rs`, `src/extract.rs`, `src/merge.rs`, and `src/mark.rs` (ANI inversion, error-rate estimation, tag coverage lookup, merge guards, and winner-take-all over-stripping protection).

## Runtime Architecture

`src/main.rs` parses the top-level CLI and dispatches to the command module matching the selected subcommand:

| Subcommand | Module | Purpose |
|------------|--------|---------|
| `extract`  | `extract` | Extract 2bRAD tags from FASTA/FASTQ into `.db` (genomes) or `.sp`/`.sylsp` (samples). |
| `merge`    | `merge`   | Merge multiple `.syldb` chunk files into one database (for parallel batch builds; refuses mixed enzyme sets). |
| `sketch`   | `sketch`  | Build k-mer sketches using sylph-style sampling into `.db`/`.sp`. |
| `query`    | `contain` | Coverage-adjusted ANI querying between databases and samples. |
| `profile`  | `contain` | Species-level taxonomic profiling with abundances and ANIs. |
| `inspect`  | `inspect` | Inspect extracted `.db`/`.sp` files. |
| `view`     | `view`    | View sketched `.db`/`.sp` files (Meta2bMap sketch format). |
| `mark`     | `mark`    | Mark taxa-specific (unique) tags in a `.db` file. |

`src/lib.rs` re-exports the public modules and the CLI struct. The crate is therefore usable both as a binary and as a library, although the public API surface is currently whatever `lib.rs` exposes.

### Module Responsibilities

- **`cmdline.rs`** — All `clap` argument structs (`ExtractArgs`, `SketchArgs`, `ContainArgs`, `ProfileArgs`, `InspectArgs`, `ViewArgs`, `MarkArgs`) and the `Mode` enum. This is the single source of truth for CLI flags.
- **`extract.rs`** — Restriction-enzyme definitions (`ENZYME_DEFINITIONS`, `ENZYME_TAG_LENGTHS`), canonical tag extraction, parallel FASTA/FASTQ processing, and serialization of `SyldbEntry` / `SylspEntry`. Also contains memory-throttling helpers.
- **`sketch.rs`** — k-mer extraction, minimap2-style hashing, cuckoo-filter deduplication, and serialization of `SequencesSketch` / `GenomeSketch`. `sketch --sylph-faithful` builds sylph-equivalent sketches (dedup instead of deleting repeated genomic k-mers, min-spacing forced to 0) for benchmarking.
- **`contain.rs`** — Core containment/ANI statistics, winner-table tag reassignment, species-level aggregation from a GTDB-style taxonomy file, and TSV output generation. This is the largest module. `profile --sketch-mode` (and `query --sketch-mode`) adapt `GenomeSketch`/`SequencesSketch` files into the same engine via `genome_sketch_to_syldb_entry` / `sketch_sample_view` (enzyme tag `kmer_c{c}_k{k}`; requires `--mismatch 0`).
- **`query.rs`** — A standalone, simpler debug-style containment query implementation that operates on `.db` and `.sp` files.
- **`inspect.rs`** / **`view.rs`** — Diagnostic commands that deserialize `.db`/`.sp` files and emit YAML or TSV summaries.
- **`mark.rs`** — Reads a `.db`, determines which tags are unique to a single genome source, and writes the file back with a `tag_uniqueness` vector per entry.
- **`inference.rs`** — Mostly commented-out statistical inference code originally intended for coverage/negative-binomial estimation.
- **`avx2_seeding.rs`** — x86_64-only AVX2 vectorized k-mer hashing, gated by `#[cfg(target_arch = "x86_64")]`. Adapted from sylph.
- **`constants.rs`** — Currently defines a `Hash = u64` alias and FNV-1a byte/string hash helpers. Much of the original content is commented out.
- **`types.rs`** — Currently almost entirely commented out; legacy shared type definitions.

### Key File Formats

- `.db` / `.sp` — Extracted tag files produced by `extract`.
- `.syldb` / `.sylsp` — Serialized `Vec<SyldbEntry>` / `Vec<SylspEntry>` files (bincode).
- The `sketch` command writes `.db`/`.sp` files in a different, sketched format inspected by `view`.

## Code Organization

```
meta2bmap/
├── Cargo.toml
├── Cargo.lock
├── README.md
└── src/
    ├── main.rs          # CLI dispatch
    ├── lib.rs           # Public re-exports
    ├── cmdline.rs       # clap CLI definitions
    ├── constants.rs     # Hash alias and FNV helpers
    ├── types.rs         # (mostly disabled legacy types)
    ├── extract.rs       # 2bRAD tag extraction
    ├── sketch.rs        # k-mer sketching
    ├── contain.rs       # ANI / profiling / taxonomy aggregation
    ├── query.rs         # Simple debug query
    ├── inspect.rs       # Inspect extracted files
    ├── view.rs          # Inspect sketched files
    ├── mark.rs          # Mark unique tags
    ├── inference.rs     # (mostly disabled stats)
    └── avx2_seeding.rs  # x86_64 AVX2 k-mer hashing
```

## Code Style Guidelines

- Use `cargo fmt` and `cargo clippy` to keep code idiomatic.
- The codebase is a mix of English and Chinese inline comments. When modifying code, prefer **English** for comments, docstrings, and user-facing messages to match the `README.md` and public CLI help text.
- Many modules contain large blocks of commented-out code (especially `types.rs`, `inference.rs`, and parts of `constants.rs`). Do not treat commented code as active logic. If you remove or replace it, update or remove the related comments.
- The project borrows algorithms and hash functions from sylph and minimap2; when adding new code that reuses such logic, preserve the existing license headers and attribution.
- Error handling uses `anyhow::Result` and `.context("...")` chains.
- Parallel loops use `rayon` (`par_iter`, `into_par_iter`). Be careful with mutable shared state; use `Mutex` / `Arc` where needed.
- Large hash maps use `FxHashMap`/`FxHashSet` for performance.

## Testing Instructions

Run the full test suite:

```bash
cargo test
```

The current suite is small (21 lib + 24 bin tests as of the sketch-mode profiling change). There are no integration tests, no CI configuration, and no sample test data in the repository.

When adding new functionality, add unit tests under `#[cfg(test)] mod tests` in the relevant module and, where practical, add a small integration test using sample FASTA/FASTQ fixtures.

## Development Workflow

1. Make changes.
2. Run `cargo fmt`.
3. Run `cargo clippy` and address warnings.
4. Run `cargo test`.
5. Build the release binary with `cargo build --release` before manual end-to-end testing.

Example end-to-end smoke test after a build:

```bash
# Extract tags from a genome FASTA
./target/release/meta2bmap extract -g genome.fa -o ./out --out-name ref_db

# Extract tags from a sample FASTQ
./target/release/meta2bmap extract -s sample_list.txt -d ./out --out-name sample

# Query
./target/release/meta2bmap query out/ref_db.db out/sample.sp > results.tsv

# Profile
./target/release/meta2bmap profile --sample-file out/sample.sp --db-file out/ref_db.db --threads 4
```

## Security Considerations

- The tool reads arbitrary user-supplied FASTA/FASTQ files and file-list text files. Paths from list files are passed directly to `File::open` / `parse_fastx_file`. Do not run with untrusted list files or paths unless the surrounding environment is trusted.
- Output files are created in user-specified directories (`-d`, `-o`, `--log-path`, etc.).
- There is no network I/O, no privileged operations, and no sandboxing; it is a local compute tool.
- `regex::bytes::Regex` is used with fixed patterns on binary sequence data; the patterns are internal constants and are not user-controlled.

## Known Project State

- This is an early-stage / research codebase. `cargo build` and `cargo test` currently succeed but produce many dead-code and unused-field warnings.
- `types.rs` and `inference.rs` are largely commented out; `constants.rs` contains a large disabled block.
- `query.rs` is a standalone debug-style implementation and is not the path used by the `query` subcommand, which is dispatched to `contain::query`.
- There is no continuous integration (`.github/workflows`, `.gitlab-ci.yml`, etc.) in the repository.
- There is no `rust-toolchain.toml` or `justfile`/`Makefile`; rely on plain Cargo commands.
