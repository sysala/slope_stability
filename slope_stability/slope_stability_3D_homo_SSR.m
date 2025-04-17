%%  Homogeneous slope and its stability  (via SSR methods)
% ======================================================================
%  This program solves a 3D slope stability problem by the modified shear
%  strength reduction method suggested in (Sysala et al. 2021). It is
%  considered the Mohr-Coulomb yield criterion, 3 Davis approaches,
%  standard finite elements (either P1 or P2 elements) and uniform meshes
%  with different densities. For P2 elements, the 7-point Gauss quadrature
%  is used. To find the safety factor of the SSR method, two continuation
%  techniques are available: the direct and the indirect techniques. The
%  indirect continuation is completed with 3 different Newton's algorithms.
%
% ======================================================================
is_agmg_present = 1;
addpath('agmg'); % for me its inside folder agmg in root of this repository  
% ======================================================================

%% The main input data
% elem_type - type of finite elements; available choices: 'P2'
elem_type = 'P2';

% Davis_type - choice of Davis' approach; available choices: 'A','B','C'
Davis_type = 'B';

%% Data from the reference element
% quadrature points and weights for volume integration
[Xi, WF] = ASSEMBLY.quadrature_volume_3D(elem_type);
% local basis functions and their derivatives
[HatP, DHatP1, DHatP2, DHatP3] = ASSEMBLY.local_basis_volume_3D(elem_type, Xi);

%% Creation/loading of the finite element mesh
% file_path = 'meshes/SSR_homo_uni.h5';
file_path = 'meshes/SSR_homo_ada_L1.h5';
% file_path = 'meshes/SSR_homo_ada_L2.h5';
% file_path = 'meshes/SSR_homo_ada_L3.h5';
% file_path = 'meshes/SSR_homo_ada_L4.h5';
% file_path = 'meshes/SSR_homo_ada_L5.h5';

switch(elem_type)
    case 'P1'
        error("Prepared meshes are only for P2 elements.")
    case 'P2'
        [coord, elem, surf, Q, material_identifier] = MESH.load_mesh_P2(file_path);
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
c0 = 6;                           % cohesion
phi = 45*pi/180;                  % frictional angle
psi = 0;                          % dilatancy angle

% elastic material parameters (FoS should be independent of these parameters)
young = 40000;                    % Young's modulus
poisson = 0.3;                    % Poisson's ratio
shear = young / (2*(1+poisson));   % shear modulus
bulk = young / (3*(1-2*poisson));  % bulk modulus
lame = bulk - 2*shear/3;           % lame's coefficient (lambda)
% specific weight of the material creating a slope
gamma = 20;

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
f_V = ASSEMBLY.vector_volume_3D(elem, coord, f_V_int, HatP, WEIGHT);

%% Input parameters for continuation (for the SSR method)
lambda_init = 0.9;        % initial lower bound of lambda
d_lambda_init = 0.1;      % initial increment of lambda
d_lambda_min = 1e-3;      % minimal increment of lambda
d_lambda_diff_scaled_min = 0.001; % minimal rate of increment of lambda
omega_max_stop = 3000;    % maximal omega, then stop
step_max = 100;           % maximal number of continuation steps

%% Input parameters for Newton's solvers
it_newt_max = 100;               % number of Newton's iterations
it_damp_max = 10;                % number of iterations within line search
tol = 1e-4;                      % relative tolerance for Newton's solvers
r_min = 1e-6;                    % basic minimal regularization of the stiffness matrix

%% Defining linear solver
linear_solver_tolerance = 1e-1;
linear_solver_maxit = 1000;
deflation_basis_tolerance = 1e-6;
linear_solver_printing = 1;

if is_agmg_present
    preconditioner_builder = @(A) LINEAR_SOLVERS.diag_prec_AGMG(A, Q);
else
    preconditioner_builder = @(A) LINEAR_SOLVERS.diag_prec_ICHOL(A, Q);
end
linear_system_solver = LINEAR_SOLVERS.DFGMRES(preconditioner_builder, ...
      linear_solver_tolerance, linear_solver_maxit, deflation_basis_tolerance, linear_solver_printing);
% linear_system_solver = LINEAR_SOLVERS.DIRECT_BACKSLASH();

%% Constitutive problem and matrix builder
dim = 3;
n_strain = dim * (dim + 1) / 2;
constitutive_matrix_builder = CONSTITUTIVE_PROBLEM.CONSTITUTIVE(B, c0, phi, psi, Davis_type, shear, bulk, lame, WEIGHT, n_strain, n_int, dim);

%--------------------------------------------------------------------------
%% Computation of the factor of safety for the SSR method

alg2on = 1; % Whether direct continuation method should be used
alg3on = 1; % Whether indirect continuation method should be used

if alg2on  % Direct continuation method - Algorithm 2
    fprintf('\n Direct continuation method - Algorithm 2\n');
    tic;
    [U2, lambda_hist2, omega_hist2, Umax_hist2] = CONTINUATION.SSR_direct_continuation(...
        lambda_init, d_lambda_init, d_lambda_min, d_lambda_diff_scaled_min, step_max, ...
        it_newt_max, it_damp_max, tol, r_min, K_elast, Q, f_V, ...
        constitutive_matrix_builder, linear_system_solver.copy());
    time_run = toc;
    fprintf("Running_time = %f \n", time_run);
end
if alg3on     % Indirect continuation method - Algorithm 3
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
    VIZ.plot_displacements_3D(U2, coord, elem);
    VIZ.plot_deviatoric_strain_3D(U2, coord, elem, B);
    % Visualization of the curve: omega -> lambda for Alg2.
    figure; hold on; box on; grid on;
    plot(omega_hist2, lambda_hist2, '-o');
    title('Direct continuation method', 'Interpreter', 'latex')
    xlabel('variable - $\xi$', 'Interpreter', 'latex');
    ylabel('strength reduction factor - $\lambda$', 'Interpreter', 'latex');
end

%% Postprocessing - visualization of selected results for ALG3
if alg3on
    VIZ.plot_displacements_3D(U3, coord, elem);
    VIZ.plot_deviatoric_strain_3D(U3, coord, elem, B);
    % Visualization of the curve: omega -> lambda for Alg3.
    figure; hold on; box on; grid on;
    plot(omega_hist3, lambda_hist3, '-o');
    title('Indirect continuation method', 'Interpreter', 'latex')
    xlabel('control variable - $\omega$', 'Interpreter', 'latex');
    ylabel('strength reduction factor - $\lambda$', 'Interpreter', 'latex');
end
