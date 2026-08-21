// src/inspect.rs

//Use this allocator when statically compiling
//instead of the default
//because the musl statically compiled binary
//uses a bad default allocator which makes the
//binary take 60% longer!!! Only affects
//static compilation though. #[cfg(target_env = "musl")]

use anyhow::Result;
use clap::Parser;
use tikv_jemallocator::Jemalloc;

mod cmdline;
mod extract;
mod sketch;
mod contain;
mod constants;
mod inspect;
mod view;
mod mark;
mod merge;

#[cfg(target_arch = "x86_64")]
mod avx2_seeding;

#[global_allocator]
static GLOBAL: Jemalloc = Jemalloc; //use std::panic::set_hook;

fn main() -> Result<()> {
    let cli = cmdline::Cli::parse();

    match cli.mode {
        cmdline::Mode::Extract(extract_args) => extract::extract(extract_args),
        cmdline::Mode::Sketch(sketch_args) => sketch::sketch(sketch_args),
        cmdline::Mode::Inspect(inspect_args) => inspect::inspect(inspect_args),
        cmdline::Mode::View(view_args) => view::view(view_args),
        cmdline::Mode::Query(contain_args) => contain::query(contain_args),
        cmdline::Mode::Profile(profile_args) => contain::profile(profile_args),
        cmdline::Mode::Mark(mark_args) => mark::mark(mark_args),
        cmdline::Mode::Merge(merge_args) => merge::merge(merge_args)
    }
}