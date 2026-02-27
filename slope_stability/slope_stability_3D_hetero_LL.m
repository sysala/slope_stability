%%  Heterogeneous slope and its stability  (via LL method)
% ======================================================================
%  This program solves a 3D slope stability problem by the limit
%  load (LL) method described in (Sysala et al., CAS 2025). The Mohr-
%  Coulomb yield criterion, Davis approach, standard finite elements 
%  (either P1 or P2 elements) and meshes with different densities are
%  considered. For P2 elements, the 11-point Gauss quadrature
%  is used. To find the safety factor of the LL method, the indirect 
%  continuation technique is used. A benchmark with a heterogeneous slope
%  from (Sysala et al., CAS 2025) is considered.
%
% ======================================================================

%% The main input data
% elem_type - type of finite elements; available choices: 'P2'
elem_type = 'P2';

% Davis_type - choice of Davis' approach; available choices: 'A','B','C'
Davis_type = 'B';

% Material parameters for each subdomain. In the following table, we
% specify in each column the following material parameters, respectively:
% [c0, phi, psi, young, poisson, gamma_sat, gamma_unsat], where
%    c0 ... Cohesion (c)
%    phi ... Friction angle (phi in degrees)
%    psi ... Dilatancy angle (psi in degrees)
%    young ... Young's modulus (E)
%    poisson ...  Poisson's ratio (nu)
%    gamma_sat ...   Specific weight - saturated (gamma_sat in kN/m^3)
%    gamma_unsat ... Specific weight - unsaturated (gamma_unsat in kN/m^3)
% If gamma_sat and gamma_unsat are not distinguished, use the same values 
% for these parameters. Each row of the table represents one subdomain. If 
% a homogeneous body is considered, only one row is prescribed.
mat_props = [15, 30,  0, 10000, 0.33, 19, 19;  % Cover layer
             15, 38,  0, 50000, 0.30, 22, 22;  % General foundation
             10, 35,  0, 50000, 0.30, 21, 21;  % Relatively weak foundation
             18, 32,  0, 20000, 0.33, 20, 20]; % General slope mass


%% Data from the reference element
% Quadrature points and weights for volume integration.
[Xi, WF] = ASSEMBLY.quadrature_volume_3D(elem_type);
% Local basis functions and their derivatives.
[HatP, DHatP1, DHatP2, DHatP3] = ASSEMBLY.local_basis_volume_3D(elem_type, Xi);


%% Creation/loading of the finite element mesh
% Available LL_hetero*.h5 meshes (nodes / elements):
%   LL_hetero_uni.h5:    23990 / 15356
%   LL_hetero_ada_L1.h5: 35867 / 22516
%   LL_hetero_ada_L2.h5: 57637 / 36504
%   LL_hetero_ada_L3.h5: 102892 / 67618
%   LL_hetero_ada_L4.h5: 191637 / 131390
%   LL_hetero_ada_L5.h5: 351693 / 246074
% file_path = 'meshes/LL_hetero_uni.h5';
file_path = 'meshes/LL_hetero_ada_L1.h5';
% file_path = 'meshes/LL_hetero_ada_L2.h5';
% file_path = 'meshes/LL_hetero_ada_L3.h5';
% file_path = 'meshes/LL_hetero_ada_L4.h5';
% file_path = 'meshes/LL_hetero_ada_L5.h5';

switch(elem_type)
    case 'P1'
        error("Prepared meshes are only for P2 elements.")
    case 'P2'
        [coord, elem, surf, Q, material_identifier] = MESH.load_mesh_P2(file_path);
        % material_identifier=ones(1,n_e) for a homogeneous body
        fprintf('P2 elements: \n')
    otherwise
        error('Bad choice of element type');
end

% Number of nodes, elements, and integration points + print.
n_n = size(coord, 2);          % Number of nodes.
n_unknown = length(coord(Q));  % Number of unknowns.
n_e = size(elem, 2);           % Number of elements.
n_q = length(WF);              % Number of quadratic points.
n_int = n_e * n_q;             % Total number of integration points.
fprintf('\n The mesh data:');
fprintf('  number of nodes = %d ', n_n);
fprintf('  number of unknowns = %d ', n_unknown);
fprintf('  number of elements = %d ', n_e);
fprintf('  number of integration points = %d \n', n_int);


%% Material parameters at integration points
% Fields with prescribed material properties
fields = {'c0',      ... % Cohesion (c)
          'phi',     ... % Friction angle (phi in degrees)
          'psi',     ... % Dilatancy angle (psi in degrees)
          'young',   ... % Young's modulus (E)
          'poisson', ... % Poisson's ratio (nu)
          'gamma_sat', ... % Specific weight - saturated (gamma_sat in kN/m^3)
          'gamma_unsat'};  % Specific weight - unsaturated (gamma_unsat in kN/m^3)

% Convert properties to structured format.
materials = cellfun(@(x) cell2struct(num2cell(x), fields, 2), num2cell(mat_props, 2), 'UniformOutput', false);

% saturation - a prescribed logical array indicating integration points 
%              where the body is saturated. If gamma_sat and gamma_unsat 
%              are the same, set saturation=true(1,n_int). Otherwise,
%              this logical array is derived from a given phreatic surface.
saturation = true(1,n_int);

% Material parameters at integration points.
[c0, phi, psi, shear, bulk, lame, gamma] = ...
      ASSEMBLY.heterogenous_materials(material_identifier, saturation, n_q, materials);


%% Assembly of the elastic matrix and volume forces vector

% Assemble the elastic stiffness matrix.
[K_elast, B, WEIGHT] = ASSEMBLY.elastic_stiffness_matrix_3D(elem, coord, shear, bulk, DHatP1, DHatP2, DHatP3, WF);

% Assemble the vector of volume forces.
% Volume forces at integration points, size (f_V_int) = (3, n_int)
f_V_int = [zeros(1, n_int); -gamma; zeros(1, n_int)];
% Compute the vector of volume forces.
f = ASSEMBLY.vector_volume_3D(elem, coord, f_V_int, HatP, WEIGHT);


%% Input parameters for the indirect continuation method

d_t_min = 1e-3;                   % Minimal increment of t.
step_max = 100;                   % Maximum number of continuation steps.
LL_omega_max = 8e7;               % Maximum value of omega.
% (Other parameters such as r_damp are not needed here.)

%% Input parameters for Newton's solvers
it_newt_max = 200;                % Number of Newton's iterations.
it_damp_max = 10;                % Number of iterations within line search.
tol = 1e-4;                      % Relative tolerance for Newton's solvers.
r_min = 1e-4;                    % Basic minimal regularization of the stiffness matrix.

%% Defining linear solver
agmg_folder = "agmg"; % Check for AGMG in specified folder
solver_type = 'DFGMRES_AGMG'; % Type of solver: "DIRECT", "DFGMRES_ICHOL", "DFGMRES_AGMG", "DFGMRES_HYPRE_BOOMERAMG"

linear_solver_tolerance = 1e-1;
linear_solver_maxit = 100;
deflation_basis_tolerance = 1e-3;
linear_solver_printing = 0;

% Optional BoomerAMG options (used when solver_type contains BOOMERAMG).
boomeramg_opts = struct('threads', 16, 'print_level', 0, ...
    'use_as_preconditioner', true);

[linear_system_solver] = LINEAR_SOLVERS.set_linear_solver(agmg_folder, solver_type, ...
    linear_solver_tolerance, linear_solver_maxit, deflation_basis_tolerance, linear_solver_printing, Q, coord, boomeramg_opts);


%% Constitutive problem and matrix builder
dim = 3;
n_strain = dim * (dim + 1) / 2;
constitutive_matrix_builder = CONSTITUTIVE_PROBLEM.CONSTITUTIVE(B, c0, phi, psi, Davis_type, shear, bulk, lame, WEIGHT, n_strain, n_int, dim);

%--------------------------------------------------------------------------
%% Computation of the limit load factor by the indirect continuation

fprintf('\n Indirect continuation method for the LL method\n');
tic;

% Compute the elastic displacement field.
constitutive_matrix_builder.reduction(1.0);
U_elast = zeros(3, n_n);              % Elastic displacement vector.
linear_system_solver.setup_preconditioner(K_elast(Q, Q));
U_elast(Q) = linear_system_solver.solve(K_elast(Q, Q), f(Q));
linear_system_solver.expand_deflation_basis(U_elast(Q));
omega_el = f(Q)' * U_elast(Q);          % Work of external forces.

% Set initial increment of omega.
d_omega_ini = omega_el / 2;
U_elast = U_elast / 2;

% Run the indirect continuation method for the LL method.
[U, t_hist, omega_hist, U_max_hist] = CONTINUATION.LL_indirect_continuation(...
    d_omega_ini, d_t_min, step_max, LL_omega_max, ...
    it_newt_max, it_damp_max, tol, r_min, K_elast, U_elast, Q, f, ...
    constitutive_matrix_builder, linear_system_solver.copy());

time_run = toc;
fprintf("Running_time = %f \n", time_run);

if contains(upper(string(solver_type)), "BOOMERAMG")
    LINEAR_SOLVERS.hypre_boomeramg_clear();
end

%% Postprocessing - visualization of selected results

VIZ.plot_displacements_3D(U, coord, elem);
VIZ.plot_deviatoric_strain_3D(U, coord, elem, B);
% Visualization of the continuation curve: omega -> t.
figure; hold on; box on; grid on;
plot(omega_hist, t_hist, '-o');
title('Indirect continuation method for the LL method', 'Interpreter', 'latex')
xlabel('control variable - $\omega$', 'Interpreter', 'latex');
ylabel('load factor - $t$', 'Interpreter', 'latex');
