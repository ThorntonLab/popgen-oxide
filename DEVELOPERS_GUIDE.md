# Developers guide

## Rust toolchain, etc...

The Rust versions used for development are specified in `rust-toolchain.toml`.
See [here](https://rust-lang.github.io/rustup/overrides.html#the-toolchain-file) for details.

We pin the version used for development because new rust versions include new lints as well as new features.
We do not want PRs to fail due to new lints, so we opt in to new style, etc., idioms by intentionally updating the version used for development.
Similarly, we don't want a PR to use a brand new language feature and break the "MSRV", which stands for "minimum supported Rust version", which is specified in `Cargo.toml`.
(Although it remains possible that a PR breaks MSRV by using features provided by the pinned toolchain version but not supported by the MSRV.
We will catch that in code review.)

## Commit message syntax

We use [conventional commits](https://www.conventionalcommits.org) for commits that should go into the change log.
In general, we only want commits in the change log that affect the primary business of the crate.
New features, code refactors, bug fixes, etc., all belong in the change log.
Changes to GitHub workflows, etc., do not need to be in the change log -- readers of the log are unlikely to care about such things.
Changes related to preparing new releases (see below) also don't go in the change log.

When working on a PR for something that is feature-gated, name the feature.
For example, for something using the `tskit` feature flag:

```
fix(tskit): Fixed some bug
```

A good way to learn our commit syntax is to read the commit log in addition to reading the URL given above.

## Version numbering

We follow [semantic versioning](https://semver.org/).
See below for how undocumented API features interact with how we treat semantic versioning.

### GitHub tags

It is arguably helpful to generate GitHub tags corresponding to releases.
For a pacakge version `x.y.z`, the corresponding tag is `vx.y.z`.
For example, `v0.31.0-alpha.2`.

NOTE: if we move to a cargo workspace, this section has to be rewritten!

## What is a breaking change?

In general, a breaking change is any change breaking the public API.
This crate uses [semver checks](https://crates.io/crates/cargo-semver-checks) during CI to check for such changes.
The semver checks tool is configured to only work on documented types and functions.
Therefore, any API change to something documented is a breaking change.

As developers, we need a way to release API features without strictly committing them to semantic versioning.
The way to do this is not to document such features.
When doing this, the new API elements **must** be decorated with `#[allow(missing_docs)]` (otherwise CI will fail).
Futher, such elements *should* have a comment block that:

* Provides minimal API documentation.
  Use `//` syntax instead of `///` so that this documentation remains a comment block and is not interpreted as a docstring.
* Contains a statement in the comment block that this feature is under development and may change in future releases.

### Version numbering and breaking changes

If CI fails because a PR contains a breaking change, **do not** bump the crate version number in that PR!
Instead, submit a separate PR that **only** bumps the version, and explain why this is happening in PR comments.
Once the version bump is merged, rebase the breaking PR.

## Deprecating API features

Deprecation can occur at any time.
Well-configured projects using this crate will start to see linting failures due to using deprecated API elements and the developers of those projects will fix that.
Removal of deprecated features should be coordinated with releases that include other breaking changes.
In theory, with breaking changes being rare, this policy will result in minimal breakage for downstream projects.

## Releasing new versions

The procedure is:

1. Update the change log and commit changes. See below for more info.
2. Update package version numbers in `Cargo.toml` files if necessary.
   (Commit the changes if necessary.)
3. Verify that the package is ready for release. See below for more info.
4. Publish the crate to crates.io.
5. Create a tag for the release.
   Be sure to have `v` at the beginning of tag names!
   For example, `git tag v2.0.0-alpha.1`.
6. Push the tag.
7. Create a release on GitHub using the tag.

### Updating the change log.

We use [git cliff](https://git-cliff.org/) to update the change log.
For example:

```
git cliff --tag v0.20.0 -u -p CHANGELOG.md
```

Then, commit the changes to `CHANGELOG.md`.
This commit **does not** use conventional syntax because we do not want it in the change log.

### Verification that the package is ready to publish

Locally, execute

```
# Clean out all non-versioned files
git clean -fd
# Validate the repo
cargo release --dry-run
```

At this point, you may need to make some commits to address any reasons why validation fails.
You need to use your best judgement about whether such commits should be in the change log or not.
Sadly, the "dry run" step does not actually catch all possible packaging errors.
(This is a known issue that the rust tools teams are working on.)

Once validated, `cargo publish` will push a release to `crates.io`.
In general, this has to be done by one of the primary developers who have write permission.

## Running tests

```sh
cargo test --all-features --all-targets
```

### Running the tests that require Python

First, make sure that [uv](https://docs.astral.sh/uv/) is installed on your system.

Then, from the directory `python/integration-tests`:

```
# Install dependencies and build the Python module
uv sync
# Run the test suite
uv run python -m pytest tests
```
