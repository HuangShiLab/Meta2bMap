//! Meta2bMap library implementation.
//!
//! The CLI binary is a thin wrapper around these modules.

pub mod cmdline;
pub mod constants;
pub mod contain;
pub mod extract;
pub mod inspect;
pub mod mark;
pub mod merge;
pub mod sketch;
pub mod view;

pub use cmdline::{Cli, Mode};
pub use constants::*;

#[cfg(target_arch = "x86_64")]
pub mod avx2_seeding;
