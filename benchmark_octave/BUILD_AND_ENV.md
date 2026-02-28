# Local High-Performance Octave Stack (OpenBLAS + librsb + sparsersb + Jupyter)

This document gives a full reproducible build and verification workflow for this repo.

Target layout (all local to repo):

- sources: `.octave_all/src`
- build trees: `.octave_all/build`
- installs: `.octave_all/install`
- octave wrapper: `.octave_all/bin/octave-rsb`
- runtime env exports: `.octave_all/env.sh`
- activation helper: `activate_optimized_octave.sh`
- jupyter venv: `.venv`

## 1. Prerequisites

You need build tools and Octave dependencies from system packages.

Typical Ubuntu/Debian baseline:

```bash
sudo apt-get update
sudo apt-get install -y \
  build-essential gfortran perl m4 pkg-config cmake ninja-build \
  bison flex gawk sed tar xz-utils curl wget \
  libreadline-dev libncurses-dev libpcre2-dev libcurl4-openssl-dev \
  libz-dev libbz2-dev liblzma-dev libfftw3-dev libarpack2-dev \
  libhdf5-dev libgraphicsmagick++1-dev libgl2ps-dev \
  libsndfile1-dev libportaudio2 \
  python3 python3-venv python3-pip
```

If some optional Octave features are missing, configure will print which libraries are absent.

## 2. Scripts Provided

Build/setup scripts are in repository root.

- `clean_local_builds.sh`
  - Removes local build/install/venv artifacts.
  - Optional:
    - `PURGE_SOURCES=1` also deletes downloaded tarballs.
    - `PURGE_RESULTS=1` also deletes benchmark JSON results.
- `build_octave_stack.sh`
  - Downloads/builds:
    - OpenBLAS `0.3.31` (`TARGET=ZEN`, OpenMP, static arch)
    - Octave `11.1.0` linked to local OpenBLAS
    - librsb `1.3.0.2` (with GCC15 compatibility patch)
    - sparsersb `1.0.9` (Octave 11 compatibility patch + `-std=gnu++17`)
  - Writes:
    - `.octave_all/env.sh`
    - `.octave_all/bin/octave-rsb`
    - `activate_optimized_octave.sh`
- `setup_jupyter_octave_venv.sh`
  - Creates local venv and installs:
    - `jupyterlab`, `notebook`, `ipykernel`, `octave_kernel`
  - Creates kernel spec:
    - name: `octave-local-rsb`
    - display: `Octave (local-rsb)`
    - uses `.octave_all/bin/octave-rsb`
- `verify_stack.sh`
  - Functional checks (`pkg load sparsersb`)
  - Performance checks at `OMP_NUM_THREADS=16`:
    - line-by-line `B'*diag(d)*B` benchmark (`full` and `sym_unique`)
    - core benchmark suite (`run_benchmarks.m`, `sparse_backend='sparsersb'`)
- `bootstrap_all.sh`
  - Runs clean -> build -> venv setup -> verification end-to-end.
  - Options:
    - `--no-clean`
    - `--no-verify`

## 3. One-Command End-to-End Setup

From repo root:

```bash
./bootstrap_all.sh
```

This is the recommended path for a fresh rebuild.

## 4. Manual Step-by-Step

```bash
# 1) Clean old local artifacts
./clean_local_builds.sh

# 2) Build OpenBLAS + Octave + librsb + sparsersb
./build_octave_stack.sh

# 3) Setup local Jupyter venv and Octave kernel
./setup_jupyter_octave_venv.sh

# 4) Verify functionality + performance envelope
THREADS=16 ./verify_stack.sh
```

## 5. Runtime Environment

To use the built stack in shell sessions:

```bash
source .octave_all/env.sh
```

or via the root activation helper:

```bash
source ./activate_optimized_octave.sh
```

Then:

```bash
"$OCTAVE_BIN" --quiet --eval "pkg load sparsersb; A=sparsersb([1;2],[1;2],[1;1],2,2); disp(full(A));"
```

or directly:

```bash
.octave_all/bin/octave-rsb --quiet
```

## 6. Jupyter Usage

```bash
source .venv/bin/activate
jupyter lab
```

Choose kernel: `Octave (local-rsb)`.

The kernel is configured to use:

- local wrapper: `.octave_all/bin/octave-rsb`
- `OMP_NUM_THREADS=16`
- local `LD_LIBRARY_PATH` with librsb and OpenBLAS

## 7. Expected Performance Envelope (this machine)

For `run_btDb_prebuilt_line_breakdown.m`, profile `medium`, `OMP_NUM_THREADS=16`:

- full constructor path (`Ifull/Jfull/vals_full`): around `0.15-0.17 s` total per 12 assemblies
- symmetric unique path (`Iu/Ju/vals_u`, `"unique","sym"`): around `0.07-0.09 s` total per 12 assemblies
- expected speedup (`full/sym_unique`): around `2.0x`

`verify_stack.sh` enforces guardrails around these ranges.

## 8. Notes on Applied Compatibility Patches

The build script patches source trees during build:

1. `librsb`:
   - moves `omp.h` include out of `extern "C"` region in `rsb_common.h` for modern GCC C++ headers.
2. `sparsersb`:
   - allows overriding hardcoded `-std=gnu++11` in `src/configure`
   - replaces deprecated Octave APIs (`mex_get_ir/jc`) with `ridx/cidx`.

These patches are reapplied automatically on each clean rebuild.

## 9. Last Full Validation Run (2026-02-28)

Validated by running:

```bash
./bootstrap_all.sh --no-clean
```

Key results (`THREADS=16`):

| Benchmark | Result |
|---|---:|
| `BtDB` prebuilt line timing (`full`) total / 12 assemblies | `0.158803 s` |
| `BtDB` prebuilt line timing (`sym_unique`) total / 12 assemblies | `0.088379 s` |
| `sym_unique` speedup vs `full` | `1.797x` |
| `B' * D_p * B` (`sparsersb` backend suite) | `0.379837 s` |
| `A * x` sparse | `0.000968 s` |
| `X' * X` dense Gram | `0.004177 s` |

Verification status: `PASSED`.

Result files written by verification:

- `benchmark_octave/results/octave_btDb_prebuilt_line_breakdown_medium_verify_t16_full.json`
- `benchmark_octave/results/octave_btDb_prebuilt_line_breakdown_medium_verify_t16_sym.json`
- `benchmark_octave/results/octave_medium_verify_sparsersb_t16.json`
