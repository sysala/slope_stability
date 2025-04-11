%%  Heterogeneous slope and its stability (via SSR methods)
% ======================================================================
%  This program solves a 3D slope stability problem for a heterogeneous slope 
%  by the modified shear strength reduction method suggested in (Sysala et al. 2021). 
%  The Mohr-Coulomb yield criterion is considered together with 3 Davis approaches. 
%  Standard finite elements (only P2 elements are available in this example) are used 
%  on uniformly discretized meshes with varying densities. For P2 elements, a 7-point 
%  Gauss quadrature is employed. To determine the safety factor of the SSR method, 
%  two continuation techniques are implemented: the direct and the indirect methods. 
%  The heterogeneous nature of the slope is modeled by assigning different material 
%  properties (e.g., cohesion, friction angle, etc.) to different domains.
%
% ======================================================================
is_agmg_present = 1;
addpath('agmg'); % For me it's inside the folder "agmg" in the root of this repository.
% ======================================================================


%% The main input data
% elem_type - type of finite elements; available choices: 'P2'
elem_type = 'P2';

% Davis_type - choice of Davis' approach; available choices: 'A', 'B', 'C'
Davis_type = 'B';


%% Data from the reference element
% Quadrature points and weights for volume integration.
% [Xi, WF] = ASSEMBLY.quadrature_volume_3D(elem_type);
[Xi, WF] = ASSEMBLY.quadrature_volume_3D_degree(8);
% Local basis functions and their derivatives.
[HatP, DHatP1, DHatP2, DHatP3] = ASSEMBLY.local_basis_volume_3D(elem_type, Xi);


%% Creation/loading of the finite element mesh
% file_path = 'meshes/SSR_hetero_uni.h5';
file_path = 'meshes/SSR_hetero_ada_L1.h5';
% file_path = 'meshes/SSR_hetero_ada_L2.h5';
% file_path = 'meshes/SSR_hetero_ada_L3.h5';
% file_path = 'meshes/SSR_hetero_ada_L4.h5';
% file_path = 'meshes/SSR_hetero_ada_L5.h5';
%file_path = 'test_michalec_opraveny2.h5';


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
% (Unified treatment for heterogeneous slopes)
% Define material properties for each domain.

% Material fields: cohesion, friction angle (deg), dilatancy angle (deg),
% Young's modulus, Poisson's ratio, and unit weight.
fields = {'c0', 'phi', 'psi', 'young', 'poisson', 'gamma'};

% Material properties for different domains:
% [c0, phi, psi, young, poisson, gamma]
mat_props = [15, 30,  0, 10000, 0.33, 19;  % Cover layer
             15, 38,  0, 50000, 0.30, 22;  % General foundation
             10, 35,  0, 50000, 0.30, 21;  % Relatively weak foundation
             18, 32,  0, 20000, 0.33, 20]; % General slope mass

% Convert properties to structured format.
materials = cellfun(@(x) cell2struct(num2cell(x), fields, 2), num2cell(mat_props, 2), 'UniformOutput', false);

% Compute material parameters at integration points.
[c0, phi, psi, shear, bulk, lame, gamma] = ASSEMBLY.heterogenous_materials(material_identifier, n_q, materials);


%% Assembly of the elastic matrix and volume forces vector

% Assemble the elastic stiffness matrix.
[K_elast, B, WEIGHT] = ASSEMBLY.elastic_stiffness_matrix_3D(elem, coord, shear, bulk, DHatP1, DHatP2, DHatP3, WF);

% Assemble the vector of volume forces.
% Volume forces at integration points, size (f_V_int) = (3, n_int)
f_V_int = [zeros(1, n_int); -gamma; zeros(1, n_int)];
% Compute the vector of volume forces.
f_V = ASSEMBLY.vector_volume_3D(elem, coord, f_V_int, HatP, WEIGHT);


%% Input parameters for continuation (for the SSR method)

lambda_init = 1.0;              % Initial lower bound of lambda.
d_lambda_init = 0.1;            % Initial increment of lambda.
d_lambda_min = 1e-3;            % Minimal increment of lambda.
d_lambda_diff_scaled_min = 0.001;% Minimal rate of increment of lambda.
omega_max_stop = 1.20E+07; %9.7969e+06;    % Maximum omega, then stop.s
step_max = 100;                 % Maximum number of continuation steps.

%% Input parameters for Newton's solvers
it_newt_max = 200;               % Number of Newton's iterations.
it_damp_max = 10;               % Number of iterations within line search.
tol = 1e-6;                     % Relative tolerance for Newton's solvers.
r_min = 1e-6;                   % Basic minimal regularization of the stiffness matrix.

%% Defining linear solver
linear_solver_tolerance = 1e-1;
linear_solver_maxit = 100;
deflation_basis_tolerance = 1e-3;
linear_solver_printing = 1;

if is_agmg_present
    preconditioner_builder = @(A) LINEAR_SOLVERS.diag_prec_AGMG(A, Q);
else
    preconditioner_builder = @(A) LINEAR_SOLVERS.diag_prec_ICHOL(A, Q);
end
linear_system_solver = LINEAR_SOLVERS.DFGMRES(preconditioner_builder, ...
     linear_solver_tolerance, linear_solver_maxit, deflation_basis_tolerance, linear_solver_printing);
%linear_system_solver = LINEAR_SOLVERS.DIRECT_BACKSLASH();

%% Constitutive problem and matrix builder
dim = 3;
n_strain = dim * (dim + 1) / 2;
constitutive_matrix_builder = CONSTITUTIVE_PROBLEM.CONSTITUTIVE(B, c0, phi, psi, Davis_type, shear, bulk, lame, WEIGHT, n_strain, n_int, dim);

%--------------------------------------------------------------------------
%% Computation of the factor of safety for the SSR method

alg2on = 0; % Use direct continuation method (Algorithm 2).
alg3on = 1; % Use indirect continuation method (Algorithm 3).

if alg2on  % Direct continuation method - Algorithm 2.
    fprintf('\n Direct continuation method - Algorithm 2\n');
    tic;
    [U2, lambda_hist2, omega_hist2, Umax_hist2, work_hist] = CONTINUATION.SSR_direct_continuation(...
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

%% Postprocessing - visualization of selected results for Algorithm 2
if alg2on
    VIZ.plot_displacements_3D(U2, coord, elem);
    VIZ.plot_deviatoric_strain_3D(U2, coord, elem, B);
    % Visualization of the curve: omega -> lambda for Algorithm 2.
    figure; hold on; box on; grid on;
    plot(omega_hist2, lambda_hist2, '-o');
    title('Direct continuation method', 'Interpreter', 'latex')
    xlabel('variable - $\xi$', 'Interpreter', 'latex');
    ylabel('strength reduction factor - $\lambda$', 'Interpreter', 'latex');
end

%% Postprocessing - visualization of selected results for Algorithm 3
if alg3on
    VIZ.plot_displacements_3D(U3, coord, elem);
    VIZ.plot_deviatoric_strain_3D(U3, coord, elem, B);
    % Visualization of the curve: omega -> lambda for Algorithm 3.
    figure; hold on; box on; grid on;
    plot(omega_hist3, lambda_hist3, '-o');
    title('Indirect continuation method', 'Interpreter', 'latex')
    xlabel('control variable - $\omega$', 'Interpreter', 'latex');
    ylabel('strength reduction factor - $\lambda$', 'Interpreter', 'latex');
end
