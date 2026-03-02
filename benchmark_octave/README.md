# benchmark_octave

Benchmarks for expensive linear-algebra kernels that mirror the slope stability code paths:

- `K_tangent = B' * D_p * B` (sparse triple product with FEM-like sparsity)
- `y = A * x` (large sparse matvec)
- `G = X' * X` (dense column dot-product matrix / Gram matrix)

For a full reproducible local build + Jupyter setup workflow, see `benchmark_octave/BUILD_AND_ENV.md`.
For `sparsersb` usage and prebuilt `B' * diag(d) * B` assembly, see `benchmark_octave/SPARSERSB_PREBUILT_GUIDE.md`.
Build/setup scripts are in `setup/` (`build_octave_stack.sh`, `setup_jupyter_octave_venv.sh`, `verify_stack.sh`), with `bootstrap_all.sh` in the repository root.

## Files

- `run_benchmarks.m`: runs one engine (MATLAB or Octave) and saves JSON results.
- `run_octave_vs_matlab.sh`: runs both engines and prints a direct comparison table.
- `compare_results.py`: compares saved JSON results.

## Usage

From repo root:

```bash
./benchmark_octave/run_octave_vs_matlab.sh medium 3 1
```

Arguments:

1. profile: `small`, `medium`, `large`
2. repeats: number of timed runs per kernel (default `3`)
3. warmup: warmup runs per kernel (default `1`)

Results are written to:

- `benchmark_octave/results/octave_<profile>.json`
- `benchmark_octave/results/matlab_<profile>.json`

## Local Octave Build

If system Octave is not installed, the runner can use a local build:

- default local path searched by script:
  - `.octave_all/install/octave-11.1.0-zen/bin/octave-cli`
- optional overrides:
  - `OCTAVE_BIN=/absolute/path/to/octave-cli`
  - `MATLAB_BIN=/absolute/path/to/matlab`
  - `OCTAVE_OMP_NUM_THREADS=<n>`

Example:

```bash
OCTAVE_OMP_NUM_THREADS=16 ./benchmark_octave/run_octave_vs_matlab.sh medium 3 1
```

Activate the local optimized build in your current shell:

```bash
source setup/activate_optimized_octave.sh
```

## Build Notes (2026-02-28)

Performed local source builds under `.octave_all/`:

1. `OpenBLAS 0.3.31` (`TARGET=ZEN`, `DYNAMIC_ARCH=0`, `USE_OPENMP=1`, `NO_AFFINITY=1`)
2. `Octave 11.1.0` built against local OpenBLAS with:
   - `CFLAGS/CXXFLAGS/FFLAGS/FCFLAGS = -O3 -march=native -mtune=native`
   - `--without-qt --without-opengl --without-fltk --disable-docs`
   - explicit `--with-blas` and `--with-lapack` to local OpenBLAS

Also tested alternate BLAS target:

- `OpenBLAS 0.3.31` (`TARGET=COOPERLAKE`) used via `LD_LIBRARY_PATH` without rebuilding Octave.

## Benchmark History

### Baseline (before rebuild)

- timestamp: `2026-02-28 14:03` (from `octave_medium.json` / `matlab_medium.json`)
- profile: `medium`, repeats `3`, warmup `1`

| Kernel                | Octave (s) | MATLAB (s) | Oct/MAT |
| --------------------- | ---------: | ---------: | ------: |
| `B' * D_p * B`        |   0.664624 |   0.054647 |  12.162 |
| `A * x (sparse)`      |   0.006833 |   0.001925 |   3.550 |
| `X' * X (dense Gram)` |   0.007076 |   0.006003 |   1.179 |

### Rebuilt Octave (local OpenBLAS ZEN, tuned `OMP_NUM_THREADS=16`)

- Octave result file: `benchmark_octave/results/octave_medium_zen_t16_seq.json`
- MATLAB comparison file: `benchmark_octave/results/matlab_medium_retest3.json`
- timestamps: `2026-02-28 14:34` (Octave), `2026-02-28 14:35` (MATLAB)

| Kernel                | Octave (s) | MATLAB (s) | Oct/MAT |
| --------------------- | ---------: | ---------: | ------: |
| `B' * D_p * B`        |   0.701838 |   0.065579 |  10.702 |
| `A * x (sparse)`      |   0.007012 |   0.002883 |   2.432 |
| `X' * X (dense Gram)` |   0.004170 |   0.008515 |   0.490 |

### Change vs baseline (Octave only)

- `B' * D_p * B`: `+5.60%` slower
- `A * x (sparse)`: `+2.62%` slower
- `X' * X (dense Gram)`: `-41.07%` faster

Interpretation:

- Dense Gram (`X' * X`) improved strongly.
- Main sparse tangent-like hotspot (`B' * D_p * B`) did not improve with this rebuild and was slightly slower.

## librsb + sparsersb Notes (2026-02-28)

Built and tested `librsb` and Octave package `sparsersb` with local Octave `11.1.0`.

### Built artifacts

- `librsb` source: `.octave_all/src/librsb-1.3.0.2.tar.gz`
- `librsb` install prefix: `.octave_all/install/librsb-1.3.0.2`
- `sparsersb` source (upstream): `.octave_all/src/sparsersb-1.0.9.tar.gz`
- `sparsersb` source (patched for Octave 11 / GCC 15): `.octave_all/src/sparsersb-1.0.9-oct11.tar.gz`

### Compatibility patches applied

1. `librsb-1.3.0.2` + GCC 15:
   - `omp.h` C++ templates conflicted with `extern "C"` in `rsb_common.h`.
   - Patched local source at `.octave_all/src/librsb-1.3.0.2/rsb_common.h` to include `omp.h` outside C linkage in C++ builds.
2. `sparsersb-1.0.9` + Octave 11:
   - package hardcoded `-std=gnu++11` (too old for Octave 11 headers).
   - patched `src/configure` to respect externally provided `SPARSERSB_CXX11`, then used `-std=gnu++17`.
   - patched `src/sparsersb.cc` to replace `SparseBoolMatrix::mex_get_ir/jc` with `ridx/cidx`.

### Install/test command shape

```bash
export PATH=.octave_all/install/librsb-1.3.0.2/bin:$PATH
export LD_LIBRARY_PATH=.octave_all/install/librsb-1.3.0.2/lib:$LD_LIBRARY_PATH
export SPARSERSB_CXX11=-std=gnu++17
.octave_all/install/octave-11.1.0-zen/bin/octave-cli --eval "pkg install -verbose .octave_all/src/sparsersb-1.0.9-oct11.tar.gz"
```

Note: `pkg load sparsersb` needs `librsb.so.0` at runtime. Keep `LD_LIBRARY_PATH` set (or install `librsb` in a system library path).

### Quick runtime checks

Snippet requested by user works:

```octave
pkg load sparsersb
A_rsb = sparsersb(A);
B_rsb = sparsersb(B);
y = A_rsb * x;
C = A_rsb * B_rsb;
```

Comparison to MATLAB using `sparsersb` backend (`medium`, repeats `3`, warmup `1`, `OMP_NUM_THREADS=16`):

| Kernel                | `sparsersb` (s) | MATLAB (s) | RSB/MAT |
| --------------------- | --------------: | ---------: | ------: |
| `B' * D_p * B`        |        0.367581 |   0.065579 |   5.605 |
| `A * x (sparse)`      |        0.001019 |   0.002883 |   0.353 |
| `X' * X (dense Gram)` |        0.004201 |   0.008515 |   0.493 |

Result files: `benchmark_octave/results/octave_medium_sparsersb_t16.json` vs `benchmark_octave/results/matlab_medium_retest3.json`.

Observed timings on this machine (`OMP_NUM_THREADS=16`):

| Kernel                                        | Octave sparse (s) | `sparsersb` (s) | RSB/Octave | Relative error |
| --------------------------------------------- | ----------------: | --------------: | ---------: | -------------: |
| `A*x` SpMV (`n=100000`, `nnz=1,000,000`)      |          0.002619 |        0.000774 |      0.296 |      4.964e-17 |
| `A*B` SpMM (`n=10000`, `nnz(A)=nnz(B)=60000`) |          0.009239 |        0.018005 |      1.949 |      0.000e+00 |

Takeaway: on this setup, `sparsersb` improves large SpMV, while SpMM was slower in the tested random case.

## COMSOL Matrix-Free vs Explicit Benchmark (2026-02-28)

Benchmark script: `benchmark_octave/run_comsol_matfree_benchmark.m`

Setup used:

- mesh path from `slope_stability_3D_hetero_seepage_SSR_comsol.m`: `slope_stability/meshes/comsol_mesh.h5`
- `element_limit = 8000` (subset of the COMSOL mesh elements)
- `matvec_repeats = 30`, `repeats = 3`, `warmup = 1`
- matrices built from the same 3D assembly path (`B`, constitutive `D_p`)

Result files:

- MATLAB: `benchmark_octave/results/matlab_comsol_matfree_e8000.json`
- Octave: `benchmark_octave/results/octave_comsol_matfree_e8000.json`
- MATLAB retest: `benchmark_octave/results/matlab_comsol_matfree_e8000_retest.json`

Requested comparison (MATLAB vs Octave, explicit path):

| Kernel                 | Octave (s) | MATLAB (s) | Oct/MAT |
| ---------------------- | ---------: | ---------: | ------: |
| `build K + 30*(K*x)`   |   1.589434 |   0.173568 |   9.157 |
| `30*(K*x), K prebuilt` |   0.138120 |   0.036354 |   3.799 |

Requested Octave mat-free run:

| Kernel                | Octave (s) |
| --------------------- | ---------: |
| `30*(B'*(D_p*(B*x)))` |   0.539164 |

Octave mat-free interpretation on this test:

- vs full explicit iteration (`build K + 30*(K*x)`): `0.339x` time (about `2.95x` faster)
- vs explicit cached-`K` matvec only (`30*(K*x)`): `3.904x` slower

Note: MATLAB printed an allocator/free error on process exit in this environment after writing JSON. A retest run produced close medians (within about 3%).

## sparsersb Input-Path Benchmark (2026-02-28)

Benchmark script: `benchmark_octave/run_sparsersb_inputpath_benchmark.m`

Run used:

```bash
OMP_NUM_THREADS=16 .octave_all/install/octave-11.1.0-zen/bin/octave-cli --quiet --eval "addpath('benchmark_octave'); run_sparsersb_inputpath_benchmark('profile','medium','repeats',3,'warmup',1,'matvec_repeats',30,'out_file','benchmark_octave/results/octave_sparsersb_inputpath_medium.json');"
```

Matrix (`medium` profile): `400000 x 400000`, `nnz = 4398920`

| Kernel                                    | Octave (s) |
| ----------------------------------------- | ---------: |
| `30*(A*x), sparse prebuilt`               |   0.201952 |
| `build sparse + 30*(A*x)`                 |   0.228747 |
| `30*(A_rsb*x), rsb prebuilt from sparse`  |   0.024118 |
| `build rsb from sparse + 30*(A_rsb*x)`    |   0.073695 |
| `build sparse + build rsb + 30*(A_rsb*x)` |   0.092564 |
| `30*(A_rsb*x), rsb prebuilt from coo`     |   0.020121 |
| `build rsb from coo + 30*(A_rsb*x)`       |   0.084470 |

Fastest paths on this setup:

- Prebuilt-call speed: `rsb prebuilt from coo` (`0.020121 s`, about `10.0x` faster than prebuilt native sparse call).
- End-to-end (include conversion + 30 calls): `build rsb from sparse + 30*(A_rsb*x)` (`0.073695 s`, about `3.10x` faster than `build sparse + 30*(A*x)`).

## Prebuilt `B'*diag(d)*B` Benchmark (2026-02-28)

Benchmark script: `benchmark_octave/run_btDb_prebuild_benchmark.m`

Goal:

- compare naive repeated assembly `A(d)=B'*diag(d)*B`
- vs optimized prebuilt-index assembly (pair pattern + `accumarray`)
- include MATLAB and Octave, and Octave `sparsersb` backend

Run configuration:

- profile: `medium`
- synthetic FEM-like `B`: `20000 x 7000`, `nnz=240000`, row nnz = `12`
- `n_d=12` (12 different `d` vectors per timed call)
- repeats `3`, warmup `1`

Result files:

- `benchmark_octave/results/matlab_btDb_prebuild_medium_native.json`
- `benchmark_octave/results/octave_btDb_prebuild_medium_native.json`
- `benchmark_octave/results/octave_btDb_prebuild_medium_sparsersb.json`

| Engine / backend       | Naive `B'*diag(d)*B` (s) | Prebuilt+`accumarray` (s) | Naive/Prebuilt |
| ---------------------- | -----------------------: | ------------------------: | -------------: |
| MATLAB / native sparse |                 0.061923 |                  0.117327 |          0.528 |
| Octave / native sparse |                 0.476603 |                  0.120283 |          3.962 |
| Octave / `sparsersb`   |                 0.428576 |                  0.153819 |          2.786 |

Precompute one-time cost (pattern build from `B`):

- MATLAB native: `0.165925 s`
- Octave native: `2.973305 s`
- Octave `sparsersb`: `3.011685 s`

Interpretation:

- In this medium setup, the prebuilt strategy is a strong win in Octave (native and `sparsersb`).
- In MATLAB, built-in sparse multiplication is faster than this prebuilt method for this case.

## Prebuilt Line-by-Line Timing (`sparsersb`) (2026-02-28)

Benchmark script: `benchmark_octave/run_btDb_prebuilt_line_breakdown.m`

Timed code block per `d`:

```octave
vals_pair = P.Wraw .* d(P.Rraw);
vals_u = accumarray(P.g, vals_pair, [numel(P.Iu), 1]);
vals_full = [vals_u; vals_u(P.off)];
A = sparsersb(P.Ifull, P.Jfull, vals_full, P.n, P.n);
```

Run configuration:

- engine: local Octave `11.1.0` + `sparsersb`
- threads: `OMP_NUM_THREADS=16`
- profile: `medium`
- constructor mode: `full`
- repeats: `10`, warmup: `2`, `n_d=12`
- matrix `B`: `20000 x 7000`, `nnz=240000`
- precompute cost (`P` build): `2.990681 s`

Result file:

- `benchmark_octave/results/octave_btDb_prebuilt_line_breakdown_medium_r10_t16_full_seq.json`

Median timing breakdown (per repeat = 12 `d` vectors):

| Line                                                    |   Median (s) |  Share (%) |
| ------------------------------------------------------- | -----------: | ---------: |
| `vals_pair = P.Wraw .* d(P.Rraw)`                       |     0.030048 |      19.16 |
| `vals_u = accumarray(P.g, vals_pair, [numel(P.Iu), 1])` |     0.013838 |       8.82 |
| `vals_full = [vals_u; vals_u(P.off)]`                   |     0.005756 |       3.67 |
| `A = sparsersb(P.Ifull, P.Jfull, vals_full, P.n, P.n)`  |     0.107183 |      68.35 |
| **TOTAL**                                               | **0.156825** | **100.00** |

Median per single `d`:

| Line                                                    | Median per `d` (s) |
| ------------------------------------------------------- | -----------------: |
| `vals_pair = P.Wraw .* d(P.Rraw)`                       |           0.002504 |
| `vals_u = accumarray(P.g, vals_pair, [numel(P.Iu), 1])` |           0.001153 |
| `vals_full = [vals_u; vals_u(P.off)]`                   |           0.000480 |
| `A = sparsersb(P.Ifull, P.Jfull, vals_full, P.n, P.n)`  |           0.008932 |

### Symmetric Constructor Hint Comparison (OMP16)

New timed constructor path:

```octave
A = sparsersb(P.Iu, P.Ju, vals_u, P.n, P.n, "unique", "sym");
```

Run configuration:

- same as above, but constructor mode: `sym_unique`
- result file: `benchmark_octave/results/octave_btDb_prebuilt_line_breakdown_medium_r10_t16_sym_seq.json`
- precompute cost (`P` build): `2.974764 s`

Median timing breakdown (per repeat = 12 `d` vectors):

| Line                                                      |   Median (s) |  Share (%) |
| --------------------------------------------------------- | -----------: | ---------: |
| `vals_pair = P.Wraw .* d(P.Rraw)`                         |     0.028647 |      36.65 |
| `vals_u = accumarray(P.g, vals_pair, [numel(P.Iu), 1])`   |     0.013183 |      16.87 |
| `A = sparsersb(P.Iu, P.Ju, vals_u, ..., "unique", "sym")` |     0.036326 |      46.48 |
| **TOTAL**                                                 | **0.078156** | **100.00** |

Before vs now (`OMP_NUM_THREADS=16`, median):

| Variant                                   | Total for 12 assemblies (s) | Per assembly (s) | Speedup vs full |
| ----------------------------------------- | --------------------------: | ---------------: | --------------: |
| Full COO mirror (`Ifull/Jfull/vals_full`) |                    0.156825 |         0.013069 |          1.000x |
| Upper-tri + `"unique","sym"`              |                    0.078156 |         0.006513 |          2.006x |

Constructor-only comparison:

| Constructor line                                      | Median for 12 assemblies (s) | Per assembly (s) | Speedup |
| ----------------------------------------------------- | ---------------------------: | ---------------: | ------: |
| `sparsersb(P.Ifull, P.Jfull, vals_full, ...)`         |                     0.107183 |         0.008932 |  1.000x |
| `sparsersb(P.Iu, P.Ju, vals_u, ..., "unique", "sym")` |                     0.036326 |         0.003027 |  2.951x |

Sanity check (same generated `B,d`, random `x`):

- relative matvec difference `norm(A_full*x - A_sym*x)/norm(A_full*x) = 1.729e-16`
