# Meta2bMap — Agent Guide

Meta2bMap is a Rust CLI for coverage-adjusted ANI queries and taxonomic profiling of
native 2bRAD and shotgun metagenomic data. It supports restriction-enzyme tag extraction,
multi-enzyme databases, a sylph-style k-mer sketch mode for controlled benchmarking, and
ANI/profiling inference.

## Commands

```bash
cargo fmt --check
cargo clippy --locked --all-targets -- -D warnings
cargo test --locked
cargo build --release --locked
```

CI runs linting, tests, and release builds on Linux and macOS.

## Architecture

`src/lib.rs` is the single library root. `src/main.rs` is a thin allocator/CLI dispatcher.
The public library exposes `Cli`, `Mode`, and the implementation modules.

| Subcommand | Implementation | Purpose |
|---|---|---|
| `extract` | `extract` | Extract canonical 2bRAD tags from FASTA/FASTQ into `.syldb`/`.sylsp`. |
| `sketch` | `sketch` | Build sampled k-mer sketches for the benchmarking/compatibility mode. |
| `query` | `contain` | Report genome/sample containment and coverage-adjusted ANI. |
| `profile` | `contain` | Winner-take-all profiling, abundance, ANI, and optional species aggregation. |
| `inspect` | `inspect` | Summarize extracted tag files. |
| `view` | `view` | Summarize k-mer sketch files. |
| `mark` | `mark` | Mark taxa-specific tags. |
| `merge` | `merge` | Merge compatible chunked `.syldb` files. |

Key modules:

- `extract.rs`: enzyme definitions, scanning, canonicalization, file formats, parallel I/O.
- `sketch.rs`: scalar and AVX2 k-mer extraction, sampling, sketch serialization.
- `contain.rs`: containment, coverage model, mismatch handling, winner reassignment,
  filtering, abundance matrices, and output.
- `avx2_seeding.rs`: x86_64 SIMD path adapted from sylph; scalar code is the reference path.

## Conventions

- Keep comments, errors, and documentation in English.
- Preserve sylph and minimap2 license attribution when adapting algorithms.
- Do not add new public APIs unless they are needed outside the crate.
- Keep profile output ordering deterministic: sort by descending ANI with stable
  sample/genome/contig tie-breakers.
- Use `anyhow::Result` with contextual errors for user-facing failures.
- Validate numeric CLI ranges and thread counts before building thread pools.

## File-format compatibility

`.syldb` v2 has an 8-byte magic header, per-tag lengths, optional tag sequences/positions,
and self-describing enzyme metadata. `.syldb` v1 remains readable, but v1 lacks position
and sequence information. Older `.sylsp` files without enzyme metadata are not compatible
and must be regenerated. Never write a new format without a round-trip test and version
fallback logic.

## Testing

Unit tests cover extraction, database round trips, gzip members, mismatch lookup, the
coverage/ANI model, winner reassignment, merge behavior, and SIMD/scalar equivalence.
`tests/cli.rs` adds end-to-end native 2bRAD extraction/profiling and invalid-parameter
checks. Add tests for bug fixes before changing the implementation.
