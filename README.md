# Slope Stability Using LL and SSR Methods

This repository contains MATLAB/Octave scripts for analyzing **slope stability problems** in 2D and 3D using the **Limit Load (LL)** and **Shear Strength Reduction (SSR)** methods. The code supports various configurations of **finite elements**, **material models**, **continuation methods**, and **linear solvers**.

---

## 📁 Repository Structure

```
README.md                   ← This file
LICENSE
bootstrap_all.sh            ← Single entry-point: builds everything
setup/                      ← Build & environment scripts
  common.sh                   Shared variables and helpers
  build_octave_stack.sh       Builds Octave + OpenBLAS + librsb + sparsersb
  setup_hypre_mex.sh          Clones/builds HYPRE (OpenMP, no MPI)
  setup_jupyter_octave_venv.sh  Sets up Jupyter + Octave kernel
  clean_local_builds.sh       Removes local build artifacts
  verify_stack.sh             Runs verification benchmarks
  activate_optimized_octave.sh  (Generated) Source to activate the environment
slope_stability/            ← Main MATLAB/Octave code
  slope_stability_*.m         Demo/driver scripts (2D/3D, LL/SSR)
  SIOPT_*.m                   Scripts for SIOPT paper experiments
  *.ipynb                     Jupyter notebook demos
  +ASSEMBLY/                  Stiffness matrix assembly (+ mex/ sources)
  +CONSTITUTIVE_PROBLEM/      Material model & yield criteria (+ mex/ sources)
  +CONTINUATION/              Continuation strategies
  +LINEAR_SOLVERS/            Iterative solvers & preconditioners (+ mex/ sources)
  +MESH/                      Mesh generation (2D) and loading (3D)
  +NEWTON/                    Newton solvers for nonlinear systems
  +SEEPAGE/                   Seepage/flow subproblem
  +VIZ/                       Visualization
  agmg/                       (User-provided) AGMG library
  docs/                       Documentation (HYPRE binding, element assembly, etc.)
  meshes/                     Prepared 3D meshes in HDF5 format
  scripts/                    Auxiliary MATLAB scripts
benchmark_octave/           ← Performance benchmarks
```


---

## ⚙️ Dependencies

- **MATLAB** (recommended version R2020b or newer).
- [**AGMG** (Algebraic Multigrid)](https://agmg.eu/) – *Optional but recommended* for best performance.
  - AGMG is **not open-source**, but an **academic license is available for free** upon request.
  - After obtaining AGMG, place the files (for linux its 'agmg.m, dmtlagmg.mexa64, zmtlagmg.mexa64') into `agmg` folder in the  'slope_stability' directory of the repository.

---

## 🚀 How to Run

### 1. Choose a script

Main scripts are located directly in the 'slope_stability' root:
- `slope_stability_2D_homo_LL.m`
- `slope_stability_2D_homo_SSR.m`
- `slope_stability_2D_Kozinec_LL.m`
- `slope_stability_3D_homo_SSR.m`
- etc.

Each script corresponds to a different setting (2D/3D, homogeneous/heterogeneous, LL/SSR method).

---

### 2. AGMG Usage

AGMG usage is controlled by solver type `DFGMRES_AGMG`:

```matlab
%% Defining linear solver
agmg_folder = 'agmg'; % Make sure this folder exists and contains AGMG
solver_type = 'DFGMRES_AGMG'; % Type of solver: "DIRECT", "DFGMRES_ICHOL", "DFGMRES_AGMG"
```

If AGMG is not available set different solver, or the code falls back to `DFGMRES_ICHOL`.

---

### 3. Selecting Mesh

#### For 2D Simulations
Meshes are generated internally based on geometry parameters (`x1`, `x2`, `y1`, `y2`, `h`, etc.) using:

```matlab
[coord, elem, ...] = MESH.mesh_P1_2D(...);
% or
[coord, elem, ...] = MESH.mesh_P2_2D(...);
```

#### For 3D Simulations
Use prepared HDF5 meshes:

```matlab
file_path = 'meshes/SSR_homo_uni.h5';
[coord, elem, ...] = MESH.load_mesh_P2(file_path);
```

You can select different levels of mesh adaptivity by commenting/uncommenting the desired file path.



---
## 🔬 HYPRE BoomerAMG Usage (OpenMP)

HYPRE BoomerAMG is available in the same slope-stability workflow as AGMG through `LINEAR_SOLVERS.set_linear_solver`.

Build once (included in `bootstrap_all.sh`, or standalone):

```bash
bash setup/setup_hypre_mex.sh
```

Then choose BoomerAMG in a slope-stability script:

```matlab
%% Defining linear solver
agmg_folder = 'agmg';
solver_type = 'DFGMRES_HYPRE_BOOMERAMG'; % also accepted: "DFGMRES_BOOMERAMG", "BOOMERAMG"

linear_solver_tolerance = 1e-1;
linear_solver_maxit = 100;
deflation_basis_tolerance = 1e-3;
linear_solver_printing = 0;

boomeramg_opts = struct('threads', 16, 'print_level', 0, ...
    'use_as_preconditioner', true);

linear_system_solver = LINEAR_SOLVERS.set_linear_solver(agmg_folder, solver_type, ...
    linear_solver_tolerance, linear_solver_maxit, ...
    deflation_basis_tolerance, linear_solver_printing, Q, coord, boomeramg_opts);
```

Notes:
- In 3D mechanics problems, passing `coord` and `Q` enables automatic construction of elasticity near-null-space vectors and DOF function mapping for BoomerAMG.
- For scalar flow-type problems (e.g., Darcy/seepage), pass `coord=[]` (or leave near-null-space empty) and set scalar options as needed.
- Full binding API and option list: `slope_stability/docs/HYPRE_MATLAB_BINDING.md`.
- Historical AGMG vs HYPRE benchmark summary: `slope_stability/docs/AGMG_vs_HYPRE_BoomerAMG.md`.


---
## � GNU Octave Support (Local Optimized Build)

The code runs on **GNU Octave** (tested with 11.1) in addition to MATLAB. A fully automated build pipeline compiles an optimized local Octave stack with all numerical dependencies.

### What gets built

| Component                  | Purpose                                               |
| -------------------------- | ----------------------------------------------------- |
| **OpenBLAS** (Zen-tuned)   | Fast BLAS/LAPACK for the local Octave                 |
| **Octave 11.1**            | CLI-only build linked to local OpenBLAS               |
| **librsb + sparsersb**     | Recursive-sparse-blocks backend for faster sparse ops |
| **HYPRE** (OpenMP, no MPI) | BoomerAMG algebraic multigrid via MEX binding         |
| **Constitutive MEX**       | OpenMP C kernels for 3D stress/tangent evaluation     |
| **Assembly MEX**           | OpenMP C kernel for element-level tangent stiffness   |
| **Jupyter kernel**         | Octave kernel in a Python venv for notebook demos     |

Everything installs under `.octave_all/` and `.venv/` inside the repository — nothing is written to system directories.

### Quick start

```bash
# Build the entire stack (Octave + HYPRE + MEX files + Jupyter):
./bootstrap_all.sh

# If HYPRE / MATLAB MEX is not needed:
./bootstrap_all.sh --skip-hypre

# Skip verification benchmarks:
./bootstrap_all.sh --no-verify
```

The build takes 15–30 minutes depending on hardware (downloads ~300 MB of sources).

### Activating the environment

After a successful build, activate the local Octave in any terminal session:

```bash
source setup/activate_optimized_octave.sh
```

This sets `PATH`, `LD_LIBRARY_PATH`, and `OMP_NUM_THREADS` (default 16). The activation is required before running any Octave command or Jupyter notebook.

### Running simulations

```bash
source setup/activate_optimized_octave.sh
cd slope_stability
octave-cli slope_stability_2D_homo_SSR.m        # 2D example
octave-cli slope_stability_3D_homo_SSR.m        # 3D example
```

### Running Jupyter notebooks

```bash
source .venv/bin/activate          # activate the Jupyter venv
jupyter lab                        # open in browser
```

Open `slope_stability/slope_stability_3D_hetero_seepage_SSR_comsol_demo.ipynb` and select the **Octave (local-rsb)** kernel.

### Individual setup scripts (in `setup/`)

| Script                         | Purpose                                           |
| ------------------------------ | ------------------------------------------------- |
| `build_octave_stack.sh`        | Builds OpenBLAS → Octave → librsb → sparsersb     |
| `setup_hypre_mex.sh`           | Clones/builds HYPRE, optionally builds MATLAB MEX |
| `setup_jupyter_octave_venv.sh` | Creates `.venv/` with Jupyter + Octave kernel     |
| `clean_local_builds.sh`        | Removes all build artifacts (keeps sources)       |
| `verify_stack.sh`              | Runs timing benchmarks to verify the build        |
| `common.sh`                    | Shared configuration (versions, paths, helpers)   |

### MEX files

C/C++ MEX sources live alongside the package they belong to:

| MEX source                                              | Package                | Built binary                              |
| ------------------------------------------------------- | ---------------------- | ----------------------------------------- |
| `+ASSEMBLY/mex/assemble_K_tangent_vals.c`               | `ASSEMBLY`             | `+ASSEMBLY/assemble_K_tangent_vals.mex`   |
| `+CONSTITUTIVE_PROBLEM/mex/constitutive_problem_3D_*.c` | `CONSTITUTIVE_PROBLEM` | `+CONSTITUTIVE_PROBLEM/*.mex`             |
| `+LINEAR_SOLVERS/mex/hypre_boomeramg_mex.cpp`           | `LINEAR_SOLVERS`       | `+LINEAR_SOLVERS/hypre_boomeramg_mex.mex` |

The MEX binaries are built automatically by `bootstrap_all.sh`. To rebuild individually:

```bash
source setup/activate_optimized_octave.sh
cd slope_stability
bash ./+ASSEMBLY/mex/build_assemble_K_tangent.sh
bash ./+CONSTITUTIVE_PROBLEM/mex/build_constitutive_3D_mex.sh
octave-cli --eval "LINEAR_SOLVERS.build_hypre_boomeramg_mex();"
```


---
## �📌 Features

These codes are focused on the solution of slope stability problems in 2D and 3D by the shear strength reduction (SSR) and limit load methods. The presented solution concept build on Davis' modifications of the methods, standard finite elements and advanced continuation techniques combined with Newton-like solvers. In particular, one can choose DAVIS A, B or C approach, P1 or P2 elements and different continutation techniques. The available folders contain various 2D and 3D problems with either homogeneous or heterogeneous geometries and will be completed stepwisely. 


---
## 📝 References

The used methods and algorithms have been developed within the following papers:

S. Sysala, M. Béreš, S. Bérešová, T. Luber, Z. Michalec: Advanced continuation and iterative methods for slope stability analysis in 3D, Computers & Structures 315, 2025, 107842, https://doi.org/10.1016/j.compstruc.2025.107842

S. Sysala, M. Béreš, S. Bérešová, J. Haslinger, J. Kružík, T. Luber: Convex Optimization Problems Inspired by Geotechnical Stability Analysis. SIAM Journal on Optimization 35(3), 1993-2016, 2025, https://doi.org/10.1137/25M1723177 

J. Karátson, S. Sysala, M. Béreš: Quasi-Newton iterative solution approaches for nonsmooth elasto-plastic problems. CAMWA 178, 61-80, 2025, [doi:10.1016/j.camwa.2024.11.022](https://doi.org/10.1016/j.camwa.2024.11.022)

S. Sysala: Advanced Continuation Methods for Limit Load and Shear Strength Reduction Methods. In P. Iványi, J. Kruis, B.H.V. Topping, (Editors), "Proceedings of the Twelfth International Conference on Engineering Computational Technology", Civil-Comp Press, Edinburgh, UK, Online volume: CCC 8, Paper 9.2, 2024, [doi:10.4203/ccc.8.9.2](http://dx.doi.org/10.4203/ccc.8.9.2)

S. Sysala, F. Tschuchnigg, E. Hrubešová, Z. Michalec: Optimization variant of the shear strength reduction method and its usage for stability of embankments with unconfined seepage. Computers and Structures 281, 2023, 107033, [doi:10.1016/j.compstruc.2023.107033](https://doi.org/10.1016/j.compstruc.2023.107033)

S. Sysala, E. Hrubešová, Z. Michalec, F. Tschuchnigg: Optimization and variational principles for the shear strength reduction method. International Journal for Numerical and Analytical Methods in Geomechanics 45, 2021, pages 2388-2407, [doi:10.1002/nag.3270](https://doi.org/10.1002/nag.3270)

J. Haslinger, S. Repin, S. Sysala: Inf-sup conditions on convex cones and applications to limit load analysis. Mathematics and Mechanics of Solids 24, 2019, pages 3331-3353, [doi:10.1142/s0218202521500330](https://doi.org/10.1177/1081286519843969)

M. Čermák, S. Sysala, J. Valdman: Efficient and flexible MATLAB implementation of 2D and 3D elastoplastic problems. Applied Mathematics and Computation 355, 2019, pages 595-614, [doi:10.1016/j.amc.2019.02.054](https://doi.org/10.1016/j.amc.2019.02.054)

S. Sysala, M. Čermák, T. Ligurský: Subdifferential-based implicit return-mapping operators in Mohr-Coulomb plasticity. ZAMM 97, 2017, pages 1502-1523, [doi:10.1002/zamm.201600215](https://doi.org/10.1002/zamm.201600215)

S. Sysala, M. Čermák, T. Koudelka, J. Kruis, J. Zeman, R. Blaheta: Subdifferential-based implicit return-mapping operators in computational plasticity. ZAMM 96, 2016, pages 1318-1338, [doi:10.1002/zamm.201500305](http://dx.doi.org/10.1002/zamm.201500305)

J. Haslinger, S. Repin, S. Sysala: Guaranteed and computable bounds of the limit load for variational problems with linear growth energy functionals. Applications of Mathematics 61, 2016, pages 527-564. [doi:10.1007/s10492-016-0146-6](http://dx.doi.org/10.1007/s10492-016-0146-6)

J. Haslinger, S. Repin, S. Sysala: A reliable incremental method of computing the limit load in deformation plasticity based on compliance: Continuous and discrete setting. Journal of Computational and Applied Mathematics 303, 2016, pages 156-170, [doi:10.1016/j.cam.2016.02.035](https://doi.org/10.1016/j.cam.2016.02.035)

S. Sysala, J. Haslinger, I. Hlaváček, M. Čermák: Discretization and numerical realization of contact problems for elastic-perfectly plastic bodies. PART I - discretization, limit analysis. ZAMM 95, 2015, pages 333-353, [doi:10.1002/zamm.201300112](https://doi.org/10.1002/zamm.201300112)

M. Čermák, J. Haslinger, T. Kozubek, S. Sysala: Discretization and numerical realization of contact problems for elastic-perfectly plastic bodies. PART II - numerical realization, limit analysis. ZAMM 95, 2015, pages 1348-1371, [doi:10.1002/zamm.201400069](https://doi.org/10.1002/zamm.201400069)
