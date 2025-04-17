%%  Heterogeneous slope and its stability  (via LL method)
% ======================================================================
%  This program solves a 3D slope stability problem for a heterogeneous slope 
%  using the modified shear strength reduction method as suggested in 
%  (Sysala et al. 2021). It considers the Mohr-Coulomb yield criterion along 
%  with 3 Davis approaches, and employs standard finite elements (here, P2 elements)
%  on uniformly discretized meshes with varying densities. For P2 elements, the 
%  7-point Gauss quadrature is used for volume integration. To determine the 
%  safety factor (limit load factor) of the SSR method, an indirect continuation 
%  technique is applied.
%
% ======================================================================

%% The main input data
% elem_type - type of finite elements; available choices: 'P2'
elem_type = 'P2';

% Davis_type - choice of Davis' approach; available choices: 'A','B','C'
Davis_type = 'B';


%% Data from the reference element
% Quadrature points and weights for volume integration.
[Xi, WF] = ASSEMBLY.quadrature_volume_3D(elem_type);
% Local basis functions and their derivatives.
[HatP, DHatP1, DHatP2, DHatP3] = ASSEMBLY.local_basis_volume_3D(elem_type, Xi);


%% Creation/loading of the finite element mesh
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
% (Unified treatment for homogeneous and heterogeneous slopes)
% Define material properties for each domain.
fields = {'c0',      ... % Cohesion (c)
          'phi',     ... % Friction angle (phi in degrees)
          'psi',     ... % Dilatancy angle (psi in degrees)
          'young',   ... % Young's modulus (E)
          'poisson', ... % Poisson's ratio (nu)
          'gamma'};      % Unit weight (gamma in kN/m^3)

% Material properties: [c0, phi, psi, young, poisson, gamma]
mat_props = [15, 30,  0, 10000, 0.33, 19;  % Cover layer
             15, 38,  0, 50000, 0.30, 22;  % General foundation
             10, 35,  0, 50000, 0.30, 21;  % Relatively weak foundation
             18, 32,  0, 20000, 0.33, 20]; % General slope mass

% Convert properties to structured format.
materials = cellfun(@(x) cell2struct(num2cell(x), fields, 2), num2cell(mat_props, 2), 'UniformOutput', false);

% Material parameters at integration points.
[c0, phi, psi, shear, bulk, lame, gamma] = ASSEMBLY.heterogenous_materials(material_identifier, n_q, materials);


%% Assembly of the elastic matrix and volume forces vector

% Assemble the elastic stiffness matrix.
[K_elast, B, WEIGHT] = ASSEMBLY.elastic_stiffness_matrix_3D(elem, coord, shear, bulk, DHatP1, DHatP2, DHatP3, WF);

% Assemble the vector of volume forces.
% Volume forces at integration points, size (f_V_int) = (3, n_int)
f_V_int = [zeros(1, n_int); -gamma; zeros(1, n_int)];
% Compute the vector of volume forces.
f = ASSEMBLY.vector_volume_3D(elem, coord, f_V_int, HatP, WEIGHT);


%% Input parameters for continuation (for the SSR method)

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
solver_type = 'DFGMRES_AGMG'; % Type of solver: "DIRECT", "DFGMRES_ICHOL", "DFGMRES_AGMG"

linear_solver_tolerance = 1e-1;
linear_solver_maxit = 100;
deflation_basis_tolerance = 1e-3;
linear_solver_printing = 0;

[linear_system_solver] = LINEAR_SOLVERS.set_linear_solver(agmg_folder, solver_type, ...
    linear_solver_tolerance, linear_solver_maxit, deflation_basis_tolerance, linear_solver_printing, Q);


%% Constitutive problem and matrix builder
dim = 3;
n_strain = dim * (dim + 1) / 2;
constitutive_matrix_builder = CONSTITUTIVE_PROBLEM.CONSTITUTIVE(B, c0, phi, psi, Davis_type, shear, bulk, lame, WEIGHT, n_strain, n_int, dim);

%--------------------------------------------------------------------------
%% Computation of the factor of safety (limit load) for the SSR method

fprintf('\n Indirect continuation method\n');
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

% Run the indirect continuation method for limit load analysis.
[U, t_hist, omega_hist, U_max_hist] = CONTINUATION.LL_indirect_continuation(...
    d_omega_ini, d_t_min, step_max, LL_omega_max, ...
    it_newt_max, it_damp_max, tol, r_min, K_elast, U_elast, Q, f, ...
    constitutive_matrix_builder, linear_system_solver.copy());

time_run = toc;
fprintf("Running_time = %f \n", time_run);

%% Postprocessing - visualization of selected results

VIZ.plot_displacements_3D(U, coord, elem);
VIZ.plot_deviatoric_strain_3D(U, coord, elem, B);
% Visualization of the curve: omega -> t.
figure; hold on; box on; grid on;
plot(omega_hist, t_hist, '-o');
title('Indirect continuation method', 'Interpreter', 'latex')
xlabel('control variable - $\omega$', 'Interpreter', 'latex');
ylabel('limit load factor - $t$', 'Interpreter', 'latex');
