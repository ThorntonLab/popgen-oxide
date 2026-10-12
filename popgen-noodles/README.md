# Minimal API for handling input data via noodles

This crate provides a simple/minimal API for processing VCF/BCF records into data types that can be understood by [popgen](https://docs.rs/popgen).

## Notes

* Projects using this crate must specify the version of both this and of [popgen](https://docs.rs/popgen) in their `Cargo.toml` files
  and let the rustc dependency resolver "work its magic".
* This crate **cannot** be used in combination with [popgen-htslib](https://docs.rs/popgen-htslib).
  Both of these crates provide different implementations of the `zlib` API.
  The result is "duplicate symbol" errors when compiling tests, binaries, etc..
  Choose one or the other.
  Currently, [rust-htslib](https://docs.rs/rust-htslib) outperforms [noodles-vcf](https://docs.rs/noodles-vcf).

