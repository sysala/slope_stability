# Slope Stability Using LL and SSR Methods

This repository contains MATLAB scripts for analyzing **slope stability problems** in 2D and 3D using the **Limit Load (LL)** and **Shear Strength Reduction (SSR)** methods. The code supports various configurations of **finite elements**, **material models**, **continuation methods**, and **linear solvers**.

---

## 📁 Repository Structure (folder 'slope_stability')

- **`slope_stability_*.m`: Main driver scripts for slope simulations in different configurations.**
- `SIOPT_*.m`: Scripts for experiments presented in the SCIOPT paper. 
- **`agmg/`: Expected to contain [AGMG](https://agmg.eu/) library (see below).**
- `+ASSEMBLY`: Assembly of stiffness matrices and force vectors.
- `+CONSTITUTIVE_PROBLEM`: Material model and yield criteria.
- `+CONTINUATION`: Implementation of continuation strategies.
- `+LINEAR_SOLVERS`: Preconditioners and iterative/direct linear solvers.
- `+MESH`: Mesh generators (2D) and mesh loaders (3D).
- `+NEWTON`: Newton solvers for nonlinear systems.
- `+VIZ`: Visualization of strains, stresses, displacements, etc.
- `meshes/`: Contains prepared 3D meshes in HDF5 format.


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

If AGMG is not available set different solver, or the code defaults to using a **direct solver**.

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
## 📌 Features

These codes are focused on the solution of slope stability problems in 2D and 3D by the shear strength reduction (SSR) and limit load methods. The presented solution concept build on Davis' modifications of the methods, standard finite elements and advanced continuation techniques combined with Newton-like solvers. In particular, one can choose DAVIS A, B or C approach, P1 or P2 elements and different continutation techniques. The available folders contain various 2D and 3D problems with either homogeneous or heterogeneous geometries and will be completed stepwisely. 


---
## 📝 References

The used methods and algorithms have been developed within the following papers:

S. Sysala, M. Béreš, S. Bérešová, T. Luber: Advanced continuation and iterative methods for slope stability analysis in 3D. Submitted in 2025.

S. Sysala, M. Béreš, S. Bérešová, J. Haslinger, J. Kružík, T. Luber: Convex Optimization Problems Inspired by Geotechnical Stability Analysis. Submitted in 2025, http://arxiv.org/abs/2312.12170.

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
