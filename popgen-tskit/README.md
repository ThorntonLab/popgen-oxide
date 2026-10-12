# API for tree sequence data input

This crate provides an API for building [popgen](https://docs.rs/popgen) data structures from [tree sequences](https://tskit.dev/learn/).
We use the [rust tskit](https://docs.rs/tskit) API as a safe wrapper for the tskit C library.
Internally, we apply incremental algorithms using edge differences to obtain count data from the tree sequence.

## Notes

* Projects using this crate must specify the version of both this and of [popgen](https://docs.rs/popgen) in their `Cargo.toml` files
  and let the rustc dependency resolver "work its magic".
