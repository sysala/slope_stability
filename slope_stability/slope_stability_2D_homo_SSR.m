%%  Homogeneous slope and its stability (via SSR methods)
% =========================================================================
%
%  This program solves a 2D slope stability problem by the modified shear
%  strength reduction method suggested in (Sysala et al. 2021). It considers
%  the Mohr-Coulomb yield criterion, three Davis approaches, standard finite
%  elements (either P1 or P2 elements), and uniform meshes with different
%  densities. For P2 elements, the 7-point Gauss quadrature is used.
%  To find the safety factor of the SSR method, two continuation techniques
%  are available: the direct and the indirect techniques.
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
h = 1/1;         % Discretization parameter

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

%% Input parameters for continuation (for the SSR method)

lambda_init = 0.9;              % Initial lower bound of lambda
d_lambda_init = 0.1;            % Initial increment of lambda
d_lambda_min = 1e-5;            % Minimal increment of lambda
d_lambda_diff_scaled_min = 0.001;% Minimal rate of increment of lambda
omega_max_stop = 7e7;           % Maximum omega, then stop
step_max = 100;                 % Maximum number of continuation steps

%% Input parameters for Newton's solvers
it_newt_max = 50;               % Number of Newton's iterations
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
%% Computation of the factor of safety for the SSR method

alg2on = 1; % Use direct continuation method (Algorithm 2).
alg3on = 1; % Use indirect continuation method (Algorithm 3).

if alg2on  % Direct continuation method - Algorithm 2.
    fprintf('\n Direct continuation method - Algorithm 2\n');
    tic;
    [U2, lambda_hist2, omega_hist2, Umax_hist2] = CONTINUATION.SSR_direct_continuation(...
        lambda_init, d_lambda_init, d_lambda_min, d_lambda_diff_scaled_min, step_max, ...
        it_newt_max, it_damp_max, tol, r_min, K_elast, Q, f_V, ...
        constitutive_matrix_builder, linear_system_solver.copy());
    time_run = toc;
    fprintf("Running_time = %f \n", time_run);
end
if alg3on     % Indirect continuation method - Algorithm 3.
    fprintf('\n Indirect continuation method - Algorithm 3 \n');
    tic;
    [U3, lambda_hist3, omega_hist3, Umax_hist3] = CONTINUATION.SSR_indirect_continuation(...
        lambda_init, d_lambda_init, d_lambda_min, d_lambda_diff_scaled_min, step_max, ...
        omega_max_stop, it_newt_max, it_damp_max, tol, r_min, K_elast, Q, f_V, ...
        constitutive_matrix_builder, linear_system_solver.copy());
    time_run = toc;
    fprintf("Running_time = %f \n", time_run);
end

%% Postprocessing - visualization of selected results for ALG2
if alg2on
    VIZ.plot_deviatoric_strain_2D(U2,coord,elem,B);
    VIZ.plot_displacements_2D(U2,coord,elem);
    % Visualization of the curve: omega -> lambda for Alg2.
    figure; hold on; box on; grid on;
    plot(omega_hist2, lambda_hist2, '-o');
    title('Direct continuation method', 'Interpreter', 'latex')
    xlabel('variable - $\xi$', 'Interpreter', 'latex');
    ylabel('strength reduction factor - $\lambda$', 'Interpreter', 'latex');
end

%% Postprocessing - visualization of selected results for ALG3
if alg3on
    VIZ.plot_deviatoric_strain_2D(U3,coord,elem,B);
    VIZ.plot_displacements_2D(U3,coord,elem);
    % Visualization of the curve: omega -> lambda for Alg3.
    figure; hold on; box on; grid on;
    plot(omega_hist3, lambda_hist3, '-o');
    title('Indirect continuation method', 'Interpreter', 'latex')
    xlabel('control variable - $\omega$', 'Interpreter', 'latex');
    ylabel('strength reduction factor - $\lambda$', 'Interpreter', 'latex');
end
