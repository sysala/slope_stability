# Element-Level Tangent Assembly (`B' * D * B` Replacement)

## Summary

The tangent stiffness matrix `K_t = B_Q' * D_t * B_Q` is assembled via a
**precomputed scatter-map and element-local quadrature** instead of two
global sparse matrix-matrix products (SpMM).  A C mex file with OpenMP
performs the numeric kernel.  The sparsity pattern of `K_t` is computed once;
only values are recomputed at each Newton iteration.

**Measured result on the COMSOL demo mesh (80 226 nodes, 55 895 elements):**

| Metric                             | Old (global `B'*D*B`) | New (element mex) |
| ---------------------------------- | --------------------: | ----------------: |
| `build_K_tangent_QQ_vals` per call |                11.4 s |           0.055 s |
| Speedup                            |                     — |          **208×** |
| K_r share of Newton time           |                79.4 % |             3.8 % |
| Newton step wall-clock             |               ~14.5 s |            ~2.8 s |
| Relative L2 error vs reference     |                     — |       5.4 × 10⁻¹⁶ |

---

## What was added

### Modified files

| File                                                      | Change                                                                                                                      |
| --------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------- |
| `+ASSEMBLY/elastic_stiffness_matrix_3D.m`                 | Returns three additional outputs `DPhi1, DPhi2, DPhi3` — the global basis-function derivatives already computed internally. |
| `+CONSTITUTIVE_PROBLEM/CONSTITUTIVE.m`                    | 12 new properties, 5 new methods (see below). `build_K_tangent_QQ_vals` now dispatches to the element path automatically.   |
| `test_notebook_export.m`                                  | Captures `DPhi*` outputs, calls `set_element_data`.                                                                         |
| `slope_stability_3D_hetero_seepage_SSR_comsol_demo.ipynb` | Same two-line change as the test script.                                                                                    |

### New files

| File                              | Purpose                                                |
| --------------------------------- | ------------------------------------------------------ |
| `mex/assemble_K_tangent_vals.c`   | C mex — OpenMP-parallel element-level assembly kernel. |
| `mex/build_assemble_K_tangent.sh` | One-liner build script (calls `mkoctfile --mex`).      |

---

## How it works

### The problem

The inner loop of the Newton solver needs the tangent stiffness matrix
restricted to free DOFs:

```
K_t(Q,Q) = B_Q' * D_t * B_Q
```

where:

* **B_Q** is the strain-displacement matrix restricted to free DOFs
  (3.7 M × 231 K, 51.5 M nonzeros, **fixed** — never changes).
* **D_t** is a block-diagonal matrix (3.7 M × 3.7 M, 614 845 dense 6×6 blocks)
  whose values come from the Mohr-Coulomb return mapping and **change every
  Newton iteration**.

Octave/MATLAB evaluate this as two consecutive SpMM operations.  Both
operations rediscover the sparsity pattern from scratch every call, even
though the pattern is identical across iterations (only values change).
This took **11.4 seconds per call** and consumed 79 % of Newton time.

### The fix: element-local assembly with precomputed scatter

The key observation is that `K_t` is really the sum of element contributions:

```
K_t = Σ_e Σ_q  B_{e,q}' · (w_q · D_{e,q}) · B_{e,q}
```

where `B_{e,q}` is a tiny 6 × 30 matrix (for 3D P2 tetrahedra with 10 nodes
× 3 DOFs) and `D_{e,q}` is a 6 × 6 constitutive tangent at quadrature point
`q` of element `e`.  The local product `B_{e,q}' D_{e,q} B_{e,q}` gives a
30 × 30 local stiffness contribution.

The approach has two phases:

#### Phase 1 — One-time precomputation

1. **Sparsity pattern**: build the elastic `K_elast(Q,Q)` once via the original
   `B_Q' * D_e * B_Q`.  Extract `(ref_I, ref_J)` — the fixed nonzero positions.
   This is the existing `init_K_r_pattern`.

2. **Scatter map**: for each element `e`, determine the global free-DOF indices
   of its 30 local DOFs.  For each local matrix entry `(a, b)`, look up which
   position `k` in `(ref_I, ref_J)` it maps to.  Store this as
   `scatter_map(a*30+b, e)` — an `int64` index into the output values vector.
   Entries mapping to Dirichlet DOFs get index 0 (skip).

   This is done by `build_scatter_map` and costs ~4 s for 56 K elements.
   Memory: `30² × n_e × 8 bytes` = **~402 MB** for the current mesh.

#### Phase 2 — Per-Newton-iteration numeric assembly

For each element `e`:

1. Zero a local 30 × 30 matrix `Ke`.
2. Loop over the element's quadrature points `q = 1 … n_q`:
   - Reconstruct the local `B_{e,q}` (6 × 30) from the stored `DPhi1, DPhi2,
     DPhi3` derivatives.
   - Read `D_{e,q}` from `DS(:, g)` (36 values → 6 × 6) and multiply by the
     quadrature weight `w_q`.
   - Accumulate `Ke += B_{e,q}' * (w_q · D_{e,q}) * B_{e,q}`.
3. Scatter: for each entry `k` of `Ke`, atomically add it to
   `V_tang[scatter_map(k, e)]`.

The output `V_tang` is a flat vector of values at the precomputed sparsity
pattern positions (`ref_I, ref_J`).  It is then used to form `K_r = r·K_elast
+ (1-r)·K_tangent` and to build a `sparsersb` matrix for the iterative solver.

### Why 208× faster

The old path did:
- Assemble a 3.7 M × 3.7 M block-diagonal `D_t` from triplets.
- Compute the intermediate `tmp = D_t * B_Q` (SpMM #1) — rediscovering
  the product's sparsity pattern.
- Compute `K_t = B_Q' * tmp` (SpMM #2) — rediscovering the pattern again.
- Extract values at known positions via `K_t(ref_idx)`.

The new path does:
- No global `D_t` construction.
- No SpMM pattern discovery.
- No intermediate sparse temporaries.
- Each element's local work is tiny (30 × 30) and memory-local.
- OpenMP parallelism over elements with atomic scatter.

---

## Data flow

```
elastic_stiffness_matrix_3D
  └─ returns: K_elast, B, WEIGHT, DPhi1, DPhi2, DPhi3
                                    │
                                    ▼
constitutive_matrix_builder.set_element_data(ELEM, DPhi1, DPhi2, DPhi3)
  └─ stores element connectivity + derivatives
                                    │
                                    ▼
constitutive_matrix_builder.init_K_r_pattern(Q)
  └─ builds K_elast(Q,Q), extracts (ref_I, ref_J, V_elast)
  └─ calls build_scatter_map()
       └─ builds scatter_map (n_local_dof² × n_e, int64)
                                    │
                at each Newton iteration:
                                    │
                                    ▼
constitutive_matrix_builder.build_K_tangent_QQ_vals()
  └─ dispatches to build_K_tangent_QQ_vals_element()
       └─ calls assemble_K_tangent_vals (C mex with OpenMP)
            inputs:  DPhi1, DPhi2, DPhi3, DS, WEIGHT, scatter_map, n_q, nnz_out
            output:  V_tang (nnz × 1 double)
```

---

## Constraints and limitations

### 1. 3D only (hardcoded `dim=3`, `n_strain=6`)

The C mex file has `#define NS 6` and `#define DIM 3`.  The B-matrix layout
(which rows correspond to ε₁₁, ε₂₂, ε₃₃, γ₁₂, γ₂₃, γ₁₃) is hardcoded
in the inner loop.

**To support 2D** (`dim=2`, `n_strain=3`), a separate mex or a runtime
dimension switch would be needed.  The pure-Octave fallback
(`build_K_tangent_QQ_vals_octave`) is also currently 3D-only.

### 2. P2 tetrahedra assumed (n_p = 10, n_local_dof = 30)

Stack arrays in the mex are allocated as `double Ke[900]`, `B_eq[180]`, etc.
For P1 elements (n_p = 4, n_local_dof = 12) or other element types, these
sizes are sufficient (they are upper bounds), but the scatter map must be
rebuilt and the B-matrix layout must match the DOF ordering convention:

```
local DOF for node i:  [(i-1)*3 + 1,  (i-1)*3 + 2,  (i-1)*3 + 3]
                     =  [u_x,          u_y,          u_z        ]
```

This matches the convention used in `elastic_stiffness_matrix_3D.m`
(see the `AUX1`/`AUX2`/`AUX3` index construction).

### 3. Scatter map memory: O(n_local_dof² × n_e)

| Mesh                      |    n_e | Map size (int64) |
| ------------------------- | -----: | ---------------: |
| Current (80 K nodes)      | 55 895 |       **402 MB** |
| 10× target (~800 K nodes) | ~560 K |        **~4 GB** |

This is well within the 128 GB budget for the target problem size.

For comparison, the **old** prebuilt row-pair index approach needed 14.6 GB
for the current mesh and ~146 GB for the 10× mesh (infeasible).

### 4. OpenMP atomic scatter

The scatter phase uses `#pragma omp atomic update` when multiple elements
contribute to the same global DOF.  This is correct but may cause contention
for nodes shared by many elements.

For the current mesh (mean row nnz of `K_t` = 81), atomics are not a
bottleneck — the 0.055 s assembly time is dominated by the element-local
`B' * D * B` contractions, not the scatter.

At very large scale or if profiling reveals scatter contention, alternatives
include:
- **Element colouring**: partition elements into independent sets so that no
  two elements in the same colour share a DOF.  Then each colour can be
  assembled without atomics.
- **Thread-local buffers**: each thread accumulates into its own `V_tang`
  copy, then reduce.  Cost: `n_threads × nnz(K_t) × 8 bytes`.

### 5. D_t must be non-symmetric capable

The Mohr-Coulomb consistent tangent with non-associated flow (sin_ψ ≠ sin_φ)
produces a **non-symmetric** D.  The element assembly handles the full 6 × 6
D and the full 30 × 30 local Ke — no symmetry is assumed.

The output `K_t` will be non-symmetric if the constitutive tangent is
non-symmetric.  This is correct and identical to the old `B_Q' * D_t * B_Q`
result.

### 6. DPhi1/2/3 must be stored (additional memory)

The global basis derivatives `DPhi1, DPhi2, DPhi3` (each `n_p × n_int`
double matrix) are stored in the `CONSTITUTIVE` object.  For P2 tetrahedra:

```
3 × (10 × 614 845) × 8 bytes = ~140 MB
```

For a 10× mesh: ~1.4 GB.  This is the same data already computed inside
`elastic_stiffness_matrix_3D` but previously discarded.

### 7. Requires compilation

The mex file must be compiled before use:

```bash
cd slope_stability/
bash mex/build_assemble_K_tangent.sh
```

This requires `mkoctfile` (Octave) or `mex` (MATLAB), plus a C compiler with
OpenMP support.  The build script is in `mex/build_assemble_K_tangent.sh`.

If the mex is not compiled, the code falls back to:
1. The pure-Octave element assembly (`build_K_tangent_QQ_vals_octave`) — correct
   but very slow due to m-code loops over elements.
2. The original global `B_Q' * D_t * B_Q` — if element data is not provided at
   all.

### 8. Fallback behaviour

The dispatch logic in `build_K_tangent_QQ_vals` is:

```
if elem_assembly_ready       →  element path (mex or Octave)
else                         →  global B_Q' * D_t * B_Q (original)
```

Element assembly becomes ready only after **both**:
- `set_element_data(ELEM, DPhi1, DPhi2, DPhi3)` has been called, **and**
- `init_K_r_pattern(Q)` has run (which triggers `build_scatter_map`).

If neither is called, the original global path is used unchanged.  This means
existing scripts that do not call `set_element_data` continue to work as
before with no code changes.

---

## How to enable in a new script

```matlab
% 1) Get DPhi outputs from the stiffness assembly
[K_elast, B, WEIGHT, DPhi1, DPhi2, DPhi3] = ASSEMBLY.elastic_stiffness_matrix_3D(
    elem, coord, shear, bulk, DHatP1, DHatP2, DHatP3, WF);

% 2) Create constitutive builder (as before)
constitutive_matrix_builder = CONSTITUTIVE_PROBLEM.CONSTITUTIVE(
    B, c0, phi, psi, Davis_type, shear, bulk, lame, WEIGHT, n_strain, n_int, 3);

% 3) Provide element data for fast assembly
constitutive_matrix_builder.set_element_data(elem, DPhi1, DPhi2, DPhi3);

% 4) Everything else is unchanged — init_K_r_pattern and the Newton solver
%    automatically use the element assembly path.
```

---

## Total memory budget (current mesh)

| Data                                   |                             Size |
| -------------------------------------- | -------------------------------: |
| B (global, kept for `build_F`)         |                           850 MB |
| B_Q (restricted to free DOFs)          |                           826 MB |
| DPhi1 + DPhi2 + DPhi3                  |                           140 MB |
| scatter_map (int64)                    |                           402 MB |
| DS (tangent moduli, changes each iter) |                           177 MB |
| ref_I, ref_J, V_elast_QQ (pattern)     |                          ~430 MB |
| **Total new overhead**                 | **~542 MB** (DPhi + scatter_map) |

For a **10× mesh** the new overhead scales to ~5.4 GB — well within 128 GB.

---

## Possible future improvements

1. **2D support**: add a 2D mex variant or runtime `dim` switch in the C code.
2. **Element colouring**: eliminate atomic scatter for even better parallel
   scaling on many-core machines.
3. **Fuse constitutive + assembly**: compute `D_{e,q}` and `B_{e,q}' D B`
   in the same compiled kernel, eliminating the global `DS` array entirely.
4. **Exploit D structure**: for elastic regions the 6×6 tangent has a known
   structure (2 parameters); hard-code the contraction to save flops.
5. **MATLAB support**: the current mex uses Octave's `mkoctfile`.  For MATLAB,
   compile with `mex -O CFLAGS='$CFLAGS -fopenmp' -lgomp assemble_K_tangent_vals.c`.
