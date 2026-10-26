# Installation

## Python wheels from a release

Each [release](https://github.com/helmutcarter/molentropy/releases) has prebuilt wheels for Python 3.10 to 3.12 on Linux (x86_64, aarch64), macOS (arm64, x86_64) and Windows (x86_64). pip can pick the right one from the release page:

```bash
pip install molentropy --find-links https://github.com/helmutcarter/molentropy/releases/expanded_assets/v0.1.0
```

or download the wheel for your Python and platform and `pip install` the file.

No Rust toolchain is needed. Building from source, below, gives the command-line program as well.

## Building from source

### Requirements

- Rust 1.94 or newer.
- For Python bindings: a Python environment with build tooling for `pyo3` (see below).

### Command line and Rust library

```bash
cargo build --release
```

The program is then at `target/release/molentropy`, or run it through cargo as `cargo run --release <path_to_parm7> <path_to_nc>`.

### Python bindings

The Python module is built with [maturin](https://www.maturin.rs/):

```bash
maturin develop --release              # into the active environment
maturin build --release --out dist     # or as a wheel
```
