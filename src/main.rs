//! Meta2bMap command-line entry point.

use anyhow::Result;
use clap::Parser;
use meta2bmap::{contain, extract, inspect, mark, merge, sketch, view, Mode};
use tikv_jemallocator::Jemalloc;

#[global_allocator]
static GLOBAL: Jemalloc = Jemalloc;

fn main() -> Result<()> {
    let cli = meta2bmap::Cli::parse();

    match cli.mode {
        Mode::Extract(args) => extract::extract(args),
        Mode::Sketch(args) => sketch::sketch(args),
        Mode::Inspect(args) => inspect::inspect(args),
        Mode::View(args) => view::view(args),
        Mode::Query(args) => contain::query(args),
        Mode::Profile(args) => contain::profile(args),
        Mode::Mark(args) => mark::mark(args),
        Mode::Merge(args) => merge::merge(args),
    }
}
