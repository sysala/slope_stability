# HYPRE MATLAB Binding API (`hypre_boomeramg_mex`)

This document describes the MATLAB API for the persistent HYPRE BoomerAMG binding implemented in this repository.

## Build and Setup

One-call setup from repository root:

```bash
./setup_hypre_mex.sh
```

This pulls/updates HYPRE sources, builds an OpenMP HYPRE install, and builds the MATLAB MEX binding.

Manual MATLAB-only build call (from `slope_stability` folder):

```matlab
LINEAR_SOLVERS.build_hypre_boomeramg_mex()
```

## MATLAB API

Two main MATLAB methods are exposed through `+LINEAR_SOLVERS`:

1. `LINEAR_SOLVERS.hypre_boomeramg_setup(A, null_space, opts, instance_id)`
2. `LINEAR_SOLVERS.hypre_boomeramg_apply(r, instance_id)`

Cleanup helper:

1. `LINEAR_SOLVERS.hypre_boomeramg_clear(instance_id)` or `LINEAR_SOLVERS.hypre_boomeramg_clear()`

### Persistence and Instances

- Setup state persists while MATLAB is running.
- Multiple instances can coexist using distinct `instance_id`.
- Calling setup with an existing `instance_id` reinitializes that instance.

## Function Signatures

### Setup

```matlab
info = LINEAR_SOLVERS.hypre_boomeramg_setup(A, null_space, opts, instance_id)
```

- `A`: real sparse `n x n` matrix.
- `null_space`: dense `n x k` near-null-space vectors, or `[]`.
- `opts`: struct of BoomerAMG options, or `[]` for defaults.
- `instance_id`: char/string or numeric scalar.

Returns `info` struct with:
- `instance_id`
- `n`
- `nnz`
- `num_functions`
- `num_interp_vectors`
- `setup_time_seconds`

### Apply

```matlab
z = LINEAR_SOLVERS.hypre_boomeramg_apply(r, instance_id)
```

- `r`: real vector of length `n`.
- Applies one BoomerAMG preconditioner action using the persistent instance.

## Default Settings (`opts = []`)

Defaults are tuned for elasticity-like systems:

- `threads = 16`
- `use_as_preconditioner = true` (forces `max_iter = 1`, `tol = 0`)
- `coarsen_type = 8`
- `agg_num_levels = 0`
- `interp_type = 17`
- `p_max_elmts = 4`
- `strong_threshold = 0.5`
- `relax_type = 8`
- `relax_sweeps = 1`
- `relax_coarse_type = 8`
- `cycle_type = 1`
- `nodal = 4`
- `nodal_diag = 1`
- `nodal_levels = 0`
- `keep_same_sign = 0`
- `interp_vec_variant = 2`
- `interp_vec_qmax = 4`
- `smooth_interp_vectors = 1`
- `interp_refine = 0`
- `num_functions = 3` only when `null_space` is provided and `mod(n,3)==0`, otherwise `1`

## `opts` Fields and Meaning

### General

- `threads`: OpenMP threads for HYPRE.
- `print_level`: BoomerAMG verbosity.
- `max_levels`: Max multigrid levels.
- `use_as_preconditioner`: If true, use one AMG cycle.
- `max_iter`: BoomerAMG iterations when not forced by preconditioner mode.
- `tol`: BoomerAMG tolerance when not forced by preconditioner mode.

### Coarsening and Interpolation

- `coarsen_type`
- `agg_num_levels`
- `interp_type`
- `p_max_elmts`
- `strong_threshold`
- `strong_threshold_r`
- `trunc_factor`
- `agg_interp_type`
- `agg_trunc_factor`
- `agg_p_max_elmts`
- `max_coarse_size`
- `min_coarse_size`

### Relaxation and Cycle

- `relax_type`
- `relax_sweeps`
- `relax_coarse_type`
- `cycle_type`

### Systems and Nodal Controls

- `num_functions`
- `dof_func` (vector length `n`, values in `[0, num_functions-1]`)
- `nodal`
- `nodal_diag`
- `nodal_levels`
- `keep_same_sign`

### Interpolation Vector Controls (for near-null-space vectors)

- `interp_vec_variant`
- `interp_vec_qmax`
- `interp_vec_abs_qtrunc`
- `smooth_interp_vectors`
- `interp_refine`

## Example A: Elasticity with Near Null Space

```matlab
% A: constrained elasticity matrix, b: constrained rhs
% coord, Q come from mesh/problem assembly
[Z_free, ~] = LINEAR_SOLVERS.near_null_space_elasticity_3D(coord, Q, true);
Z = zeros(size(A,1), size(Z_free,2));
Z(Q(:), :) = Z_free;

opts = struct( ...
    'threads', 16, ...
    'num_functions', 3, ...
    'nodal', 4, ...
    'nodal_diag', 1, ...
    'coarsen_type', 8, ...
    'interp_type', 17, ...
    'strong_threshold', 0.5, ...
    'interp_vec_variant', 2, ...
    'interp_vec_qmax', 4, ...
    'smooth_interp_vectors', 1, ...
    'use_as_preconditioner', true);

info = LINEAR_SOLVERS.hypre_boomeramg_setup(A, Z, opts, 'elas_l2');
M = @(r) LINEAR_SOLVERS.hypre_boomeramg_apply(r, 'elas_l2');
[x, flag, relres, iter] = pcg(A, b, 1e-8, 200, M, []);
LINEAR_SOLVERS.hypre_boomeramg_clear('elas_l2');
```

## Example B: Darcy-like SPD Problem (No Null Space)

```matlab
% A: SPD Darcy matrix, b: rhs
opts = struct( ...
    'threads', 16, ...
    'num_functions', 1, ...
    'nodal', 0, ...
    'coarsen_type', 10, ...
    'interp_type', 6, ...
    'strong_threshold', 0.25, ...
    'use_as_preconditioner', true);

info = LINEAR_SOLVERS.hypre_boomeramg_setup(A, [], opts, 'darcy');
M = @(r) LINEAR_SOLVERS.hypre_boomeramg_apply(r, 'darcy');
[x, flag, relres, iter] = pcg(A, b, 1e-8, 300, M, []);
LINEAR_SOLVERS.hypre_boomeramg_clear('darcy');
```

## Notes

- The binding is single-process (`MPI` disabled build), OpenMP threaded.
- For multiple right-hand sides, call `apply` repeatedly.
- Call `LINEAR_SOLVERS.hypre_boomeramg_clear()` before unloading mex or at end of scripts to release memory.
