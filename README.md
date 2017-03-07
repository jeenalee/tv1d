# `tv1d`

[![](http://meritbadge.herokuapp.com/tv1d)](https://crates.io/crates/tv1d) [![](https://docs.rs/tv1d/badge.svg)](https://docs.rs/tv1d/) [![Build Status](https://travis-ci.org/jeenalee/tv1d.png?branch=master)](https://travis-ci.org/jeenalee/tv1d/)

Total variation denoising algorithms for 1d data.

## Total Variation Denoising

[Total variation denoising](https://en.wikipedia.org/wiki/Total_variation_denoising) algorithms denoise signals through reducing the total variation. As they are effective at preserving edges while removing the noise, total variation denoising algorithms are especially suitable for signals whose model are [piecewise constant function](https://en.wikipedia.org/wiki/Step_function). In short, piecewise constant 1D data would look like a mix of flat regions and jumps between them.

## Install

Add this to your `Cargo.toml`:
```toml
[dependencies]
tv1d = "0.1.0"
```

## Example Usage

``` rust
extern crate tv1d

fn main() {
    let input = vec![13.0, 24.3, 63.41, 13.6];
    let lambda = 3.0;

    let output = tv1d::condat(&input, lambda);
}
```

## Documentation

Read the [documentation on Docs.rs](https://docs.rs/tv1d).

## License

This crate is licensed under MIT license ([`LICENSE`](./LICENSE)).

## Using Rust crate from other languages

Please check out the Rust Book's chapter ["Rust Inside Other Languages"](https://doc.rust-lang.org/1.2.0/book/rust-inside-other-languages.html).

## Contribution

See [CONTRIBUTING.md](./CONTRIBUTING.md)!
