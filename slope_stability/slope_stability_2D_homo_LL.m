%%  Homogeneous slope and its stability (via LL method)
% =========================================================================
%
%  This program solves a 2D slope stability problem by the limit
%  load (LL) method described in (Sysala et al., CAS 2025). The Mohr-
%  Coulomb yield criterion, Davis approach, standard finite elements 
%  (either P1 or P2 elements) and meshes with different densities are
%  considered. For P2 elements, the 7-point Gauss quadrature
%  is used. To find the safety factor of the LL method, the indirect 
%  continuation technique is used. A benchmark with a homogeneous slope
%  is considered. It is possible to change slope inclination and other
%  geometrical parameters.
%
% =========================================================================

%% Main input data

% elem_type - type of finite elements; available choices: 'P1', 'P2'
elem_type = 'P2';

% Davis_type - choice of Davis' approach; available choices: 'A', 'B', 'C'
Davis_type = 'B';

%% Geometrical parameters
x1 = 15;         % Length of the body in front of the slope
x3 = 15;         % Length of the body behind the slope
y1 = 10;         % Height of the body below the slope
y2 = 10;         % Height of the slope
beta = pi/4;     % Slope angle
x2 = y2/tan(beta); % Length of the slope in the x-direction

%% Mesh data
h = 1/2;         % Discretization parameter

%% Strength material parameters
c0 = 6;                           % Cohesion
phi = 45*pi/180;                  % Frictional angle
psi = 0;                          % Dilatancy angle

%% Elastic material parameters (FoS should be independent of these parameters)
young = 40000;                     % Young's modulus
poisson = 0.3;                     % Poisson's ratio
shear = young / (2 * (1 + poisson)); % Shear modulus
bulk = young / (3 * (1 - 2 * poisson)); % Bulk modulus
lame = bulk - 2 * shear / 3;        % Lame's coefficient (lambda)

% Specific weight of the material creating a slope
gamma = 20;

%% Data from the reference element

% Quadrature points and weights for volume integration
[Xi, WF] = ASSEMBLY.quadrature_volume_2D(elem_type);
% Local basis functions and their derivatives
[HatP, DHatP1, DHatP2] = ASSEMBLY.local_basis_volume_2D(elem_type, Xi);

%% Creation of the uniform finite element mesh

switch(elem_type)
    case 'P1'
        [coord, elem, ELEM_ED, EDGE_EL, Q] = MESH.mesh_P1_2D(h, x1, x2, x3, y1, y2);
        fprintf('P1 elements: \n')
    case 'P2'
        [coord, elem, ELEM_ED, EDGE_EL, Q] = MESH.mesh_P2_2D(h, x1, x2, x3, y1, y2);
        fprintf('P2 elements: \n')
    otherwise
        error('Bad choice of element type');
end

% Number of nodes, elements, and integration points + print
n_n = size(coord,2);          % Number of nodes
n_unknown = length(coord(Q)); % Number of unknowns
n_e = size(elem,2);           % Number of elements
n_ed = size(EDGE_EL,2);       % Number of edges
n_q = length(WF);             % Number of quadrature points
n_int = n_e * n_q;            % Total number of integration points

fprintf('\n The mesh data:');
fprintf('  Number of nodes = %d ', n_n);
fprintf('  Number of unknowns = %d ', n_unknown);
fprintf('  Number of elements = %d ', n_e);
fprintf('  Number of edges = %d ', n_ed);
fprintf('  Number of integration points = %d \n', n_int);

%% Material parameters at integration points
% (For a unified treatment of homogeneous and heterogeneous slopes)
c0 = c0 * ones(1, n_int);
phi = phi * ones(1, n_int);
psi = psi * ones(1, n_int);
shear = shear * ones(1, n_int);
bulk = bulk * ones(1, n_int);
lame = lame * ones(1, n_int);
gamma = gamma * ones(1, n_int);

%% Assembling of the elastic stiffness matrix
[K_elast, B, WEIGHT] = ASSEMBLY.elastic_stiffness_matrix_2D(elem, coord, DHatP1, DHatP2, WF, shear, lame);

%% Assembling of the vector of volume forces

% Volume forces at integration points, size(f_V_int) = (2, n_int)
f_V_int = [zeros(1, n_int); -gamma];
% Vector of volume forces
f_V = ASSEMBLY.vector_volume_2D(elem, coord, f_V_int, HatP, WEIGHT);

%% Input parameters for the indirect continuation method
d_t_min = 1e-3;                % Minimal increment of load factor t.
step_max = 100;                % Maximum number of continuation steps.
LL_omega_max = 2000;            % Maximum value of the control parameter omega.

%% Input parameters for Newton's solvers
it_newt_max = 100;               % Number of Newton's iterations
it_damp_max = 10;               % Number of iterations within line search
tol = 1e-4;                     % Relative tolerance for Newton's solvers
r_min = 1e-4;                   % Basic minimal regularization of the stiffness matrix

%% Defining linear solver
agmg_folder = "agmg"; % Check for AGMG in specified folder
solver_type = 'DIRECT'; % Type of solver: "DIRECT", "DFGMRES_ICHOL", "DFGMRES_AGMG"

linear_solver_tolerance = 1e-1;
linear_solver_maxit = 100;
deflation_basis_tolerance = 1e-3;
linear_solver_printing = 0;

[linear_system_solver] = LINEAR_SOLVERS.set_linear_solver(agmg_folder, solver_type, ...
    linear_solver_tolerance, linear_solver_maxit, deflation_basis_tolerance, linear_solver_printing, Q);


%% Constitutive problem and matrix builder
dim = 2;
n_strain = dim * (dim + 1) / 2;
constitutive_matrix_builder = CONSTITUTIVE_PROBLEM.CONSTITUTIVE(B, c0, phi, psi, Davis_type, shear, bulk, lame, WEIGHT, n_strain, n_int, dim);

%--------------------------------------------------------------------------
%% Computation of the limit load factor by the indirect continuation
fprintf('\n Indirect continuation method for the LL method\n');
tic;

% Compute the elastic displacement field as a starting point.
constitutive_matrix_builder.reduction(1.0);
U_elast = zeros(2, n_n);              % Elastic displacement vector.
linear_system_solver.setup_preconditioner(K_elast(Q, Q));
U_elast(Q) = linear_system_solver.solve(K_elast(Q, Q), f_V(Q));
linear_system_solver.expand_deflation_basis(U_elast(Q));
omega_el = f_V(Q)' * U_elast(Q);        % Work of external forces.

% Set the initial increment of omega.
d_omega_ini = omega_el / 5;
% Scale the elastic displacement field.
U_elast = U_elast / 5;

% Run the indirect continuation method for the LL method.
[U, t_hist, omega_hist, U_max_hist] = CONTINUATION.LL_indirect_continuation(...
    d_omega_ini, d_t_min, step_max, LL_omega_max, ...
    it_newt_max, it_damp_max, tol, r_min, K_elast, U_elast, Q, f_V, ...
    constitutive_matrix_builder, linear_system_solver.copy());

time_run = toc;
fprintf("Running_time = %f \n", time_run);


%% Postprocessing - visualization of selected results
VIZ.plot_deviatoric_strain_2D(U,coord,elem,B);
VIZ.plot_displacements_2D(U,coord,elem);
% Visualization of the continuation curve: omega -> lambda.
figure; hold on; box on; grid on;
plot(omega_hist, t_hist, '-o');
title('Indirect continuation method for the LL method', 'Interpreter', 'latex');
xlabel('Control variable - $\omega$', 'Interpreter', 'latex');
ylabel('Load factor - $t$', 'Interpreter', 'latex');
