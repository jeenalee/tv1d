//! TODO

#![deny(missing_docs)]
#![deny(missing_debug_implementations)]

#![cfg_attr(feature="clippy", feature(plugin))]

#![cfg_attr(feature="clippy", plugin(clippy))]

extern crate num;

mod denoise;
pub use denoise::*;

mod utils;
pub use utils::*;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {}
}
