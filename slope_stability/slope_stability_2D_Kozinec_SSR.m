%%  Heterogeneous slope in Doubrava-Kozinec and its stability (SSR methods)
% =========================================================================
%
%  This program solves a 2D slope stability problem by the modified shear
%  strength reduction method suggested in (Sysala et al. 2021). It is
%  considered the Mohr-Coulomb yield criterion, 3 Davis approaches,
%  standard finite elements (P1, P2 or P4 elements) and uniform meshes
%  with different densities. Gauss quadrature is used for numerical
%  integration. To find the safety factor of the SSR method, 2 continuation
%  techniques are available: the direct and the indirect techniques. The
%  indirect continuation is completed with 3 different Newton's algorithms.
%
% ======================================================================
%

%% The main input data

% elem_type - type of finite elements; available choices: 'P1', 'P2', 'P4'
elem_type='P2';

% Davis_type - choice of Davis' approach; available choices: 'A','B','C'
Davis_type='B';

%
% Data from the reference element
%
% quadrature points and weights for volume integration
[Xi, WF] = ASSEMBLY.quadrature_volume_2D(elem_type);
% local basis functions and their derivatives
[HatP,DHatP1,DHatP2] = ASSEMBLY.local_basis_volume_2D(elem_type, Xi);


%
%  Loading of the mesh imported from COMSOL
%
[coord, elem, Q, mater] = MESH.load_mesh_Kozinec(elem_type, 'meshes/Kozinec/');
% number of nodes, elements and integration points + print
n_n=size(coord,2);
n_unknown=length(coord(Q)); % number of unknowns
n_e=size(elem,2);           % number of elements
n_q=length(WF);             % number of quadratic points
n_int = n_e*n_q ;           % total number of integrations points
%
fprintf('\n');
fprintf('Mesh data:');
fprintf('  number of nodes =%d ',n_n);
fprintf('  number of unknowns =%d ',n_unknown);
fprintf('  number of elements =%d ',n_e);
fprintf('  number of integration points =%d ',n_int);
fprintf('\n');

%
% Material parameters at integration points
%
%% Material Parameters at Integration Points
% (Unified treatment for heterogeneous slopes)
% Define material properties for each domain.

% Elastic material parameters (FoS should be independent of these parameters)
young = 16000;                    % Young's modulus
poisson = 0.4;                     % Poisson's ratio
shear = young / (2 * (1 + poisson)); % Shear modulus
bulk = young / (3 * (1 - 2 * poisson)); % Bulk modulus
lame = bulk - 2 * shear / 3;      % Lame's coefficient (lambda)

% Material fields: cohesion, friction angle (rad), dilatation angle,
% saturated weight, unsaturated weight.
fields = {'c0', 'phi', 'psi', 'gamma_sat', 'gamma_unsat'};

% Material properties for different domains:
% Expanded list to include all subdomains explicitly
% [c0, phi (rad), psi, gamma_sat, gamma_unsat]
mat_props = [9,  26 * pi / 180, 0, 20.3, 20.7;  % Material #1
    2,  33 * pi / 180, 0, 19.0, 20.5;  % Material #2
    5,  27 * pi / 180, 0, 19.4, 21.4;  % Material #3
    3,  13 * pi / 180, 0, 20.0, 20.5;  % Material #4
    5,  27 * pi / 180, 0, 19.4, 21.4;  % Material #5 (same as #3)
    3,  13 * pi / 180, 0, 20.0, 20.5;  % Material #6 (same as #4)
    1,  45 * pi / 180, 0, 20.5, 20.6]; % Material #7

% Convert properties to structured format.
materials = cellfun(@(x) cell2struct(num2cell(x), fields, 2), num2cell(mat_props, 2), 'UniformOutput', false);

% Compute material parameters at integration points.

% Compute material parameters at integration points.
[c0, phi, psi, gamma_sat, gamma_unsat] = ASSEMBLY.heter_mater(mater, n_q, materials);

% Specific weight at integration points depending on a given saturation curve.
% This curve is prescribed within the function "gravity".
gamma = ASSEMBLY.gravity(gamma_sat, gamma_unsat, coord, elem, HatP);

% Homogeneous (elastic) parameters at integration points.
shear = shear * ones(1, n_int);
bulk = bulk * ones(1, n_int);
lame = lame * ones(1, n_int);



%
% Assembling of the elastic stiffness matrix
%
[K_elast,B,WEIGHT]=ASSEMBLY.elastic_stiffness_matrix_2D(elem,coord,...
    DHatP1,DHatP2,WF,shear,lame);

%
% Assembling of the vector of volume forces
%

% volume forces at integration points, size(f_V_int)=(2,n_int)
f_V_int = [zeros(1,n_int);-gamma] ;
% vector of volume forces
f_V=ASSEMBLY.vector_volume_2D(elem,coord,f_V_int,HatP,WEIGHT);


%% Input parameters for continuation (for the SSR method)

lambda_init = 0.7;              % Initial lower bound of lambda
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
