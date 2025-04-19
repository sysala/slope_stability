%%  Homogeneous slope and its stability  (via LL method)
% ======================================================================
%  This program solves a 3D slope stability problem by the limit
%  load method suggested in (Sysala et al. 2021). It is
%  considered the Mohr-Coulomb yield criterion, 3 Davis approaches,
%  standard finite elements (either P1 or P2 elements) and meshes
%  with different densities. For P2 elements, the 11-point Gauss quadrature
%  is used. To find the safety factor of the SSR method, two continuation
%  techniques are available: the direct and the indirect techniques.
%
% ======================================================================

%% The main input data
% elem_type - type of finite elements; available choices: 'P2'
elem_type = 'P2';

% Davis_type - choice of Davis' approach; available choices: 'A','B','C'
Davis_type = 'B';

lambda_ell = 1.0; % for Limit Load, factor of strength reduction (for ell plot)

%% Data from the reference element
% quadrature points and weights for volume integration
[Xi, WF] = ASSEMBLY.quadrature_volume_3D(elem_type);
% local basis functions and their derivatives
[HatP, DHatP1, DHatP2, DHatP3] = ASSEMBLY.local_basis_volume_3D(elem_type, Xi);

%% Creation/loading of the finite element mesh
file_path = 'meshes/SIOPT_L0.h5';
% file_path = 'meshes/SIOPT_L1.h5';
% file_path = 'meshes/SIOPT_L5.h5';

switch(elem_type)
    case 'P1'
        error("Prepared meshes are only for P2 elements.")
    case 'P2'
        [coord, elem, surf, Q, material_identifier] = MESH.load_mesh_P2(file_path, 1);
        fprintf('P2 elements: \n')
    otherwise
        error('bad choice of element type');
end

% number of nodes, elements and integration points + print
n_n = size(coord,2);          % number of nodes
n_unknown = length(coord(Q)); % number of unknowns
n_e = size(elem,2);           % number of elements
n_q = length(WF);             % number of quadratic points
n_int = n_e * n_q;            % total number of integration points
%
fprintf('\n The mesh data:');
fprintf('  number of nodes = %d ', n_n);
fprintf('  number of unknowns = %d ', n_unknown);
fprintf('  number of elements = %d ', n_e);
fprintf('  number of integration points = %d \n', n_int);

%% Material parameters at integration points
% (for a unified treatment of homogeneous and heterogeneous slopes)
% Define material properties for each domain

% strength material parameters
c0 = 15;                           % Cohesion [kPa]
phi = 20 * pi / 180;               % Friction angle [radians] (20 degrees)
psi = phi;                           % Dilation angle [radians] (assume no dilation)

% Elastic material parameters (not necessary for limit analysis)
young = 40000;                      % Young's modulus [kPa] (typical for soft clay)
poisson = 0.3;                     % Poisson's ratio (soft, near-saturated soil)
shear = young / (2 * (1 + poisson)); % Shear modulus [kPa]
bulk = young / (3 * (1 - 2 * poisson)); % Bulk modulus [kPa]
lame = bulk - 2 * shear / 3;       % Lame's coefficient (lambda) [kPa]

% Specific weight of the soil material
gamma = 20;                     % Unit weight [kN/m^3] (calculated from density)

c0 = c0 * ones(1, n_int);
phi = phi * ones(1, n_int);
psi = psi * ones(1, n_int);
shear = shear * ones(1, n_int);
bulk = bulk * ones(1, n_int);
lame = lame * ones(1, n_int);
gamma = gamma * ones(1, n_int);

%% Assembly of the elastic matrix and volume forces vector

% Assembling of the elastic stiffness matrix
[K_elast, B, WEIGHT] = ASSEMBLY.elastic_stiffness_matrix_3D(elem, coord, shear, bulk, DHatP1, DHatP2, DHatP3, WF);

% Assembling of the vector of volume forces
% volume forces at integration points, size(f_V_int) = (3, n_int)
f_V_int = [zeros(1, n_int); -gamma; zeros(1, n_int)];
% vector of volume forces
f = ASSEMBLY.vector_volume_3D(elem, coord, f_V_int, HatP, WEIGHT);

%% Input parameters for continuation (for the SSR method)
lambda_init = 0.9;        % initial lower bound of lambda
d_lambda_init = 0.1;      % initial increment of lambda
d_t_min = 1e-3;      % minimal increment of lambda
d_lambda_diff_scaled_min = 0.001; % minimal rate of increment of lambda
LL_omega_max = 3.6e7;    % maximal omega, then stop
step_max = 100;           % maximal number of continuation steps

%% Input parameters for Newton's solvers
it_newt_max = 100;               % number of Newton's iterations
it_damp_max = 10;                % number of iterations within line search
tol = 1e-4;                      % relative tolerance for Newton's solvers
r_min = 1e-4;                    % basic minimal regularization of the stiffness matrix

%% Defining linear solver
agmg_folder = "agmg"; % Check for AGMG in specified folder
solver_type = 'DFGMRES_AGMG'; % Type of solver: "DIRECT", "DFGMRES_ICHOL", "DFGMRES_AGMG"

linear_solver_tolerance = 1e-1;
linear_solver_maxit = 200;
deflation_basis_tolerance = 1e-3;
linear_solver_printing = 1;

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
constitutive_matrix_builder.reduction(lambda_ell);
U_elast = zeros(3, n_n);              % Elastic displacement vector.
linear_system_solver.setup_preconditioner(K_elast(Q, Q));
U_elast(Q) = linear_system_solver.solve(K_elast(Q, Q), f(Q));
linear_system_solver.expand_deflation_basis(U_elast(Q));
omega_el = f(Q)' * U_elast(Q);          % Work of external forces.

% Set initial increment of omega.
d_omega_ini = omega_el / 5;
U_elast = U_elast / 5;

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