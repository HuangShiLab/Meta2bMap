pub mod cmdline;
pub mod constants;
pub mod query;
pub mod extract;
pub mod sketch;
pub mod inspect;
pub mod contain;


pub use cmdline::Cli;
pub use constants::*;

#[cfg(target_arch = "x86_64")]
pub mod avx2_seeding;


