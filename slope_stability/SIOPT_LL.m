%%  Homogeneous slope and its stability  (via LL method)
% ======================================================================
%  This program solves a 2D slope stability problem by the limit
%  load (LL) method described in (Sysala et al., CAS 2025). The Mohr-
%  Coulomb yield criterion, Davis approach, standard finite elements 
%  (either P1 or P2 elements) and meshes with different densities are
%  considered. For P2 elements, the 11-point Gauss quadrature
%  is used. To find the safety factor of the LL method, the indirect 
%  continuation technique is used. The benchmark described in the paper 
%  (Sysala et al., SIOPT 2025) is considered.
%
% ======================================================================

%% The main input data
% elem_type - type of finite elements; available choices: 'P2'
elem_type = 'P2';

% Davis_type - choice of Davis' approach; available choices: 'A','B','C'
Davis_type = 'B';

lambda_ell = 1.0; % choose lambda_ell = 1.0 for the LL method, 
                  % to construct the function ell, which relates the LL and 
                  % SSR method, choose other values of lambda_ell

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
mat_props = [15, 20, 20, 40000, 0.3, 20, 20]; 

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
        [coord, elem, surf, Q, ~] = MESH.load_mesh_P2(file_path, 1);
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

% The array material_identifier for a homogeneous body
material_identifier = zeros(1,n_e);

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

% Assembling of the elastic stiffness matrix
[K_elast, B, WEIGHT] = ASSEMBLY.elastic_stiffness_matrix_3D(elem, coord, shear, bulk, DHatP1, DHatP2, DHatP3, WF);

% Assembling of the vector of volume forces
% volume forces at integration points, size(f_V_int) = (3, n_int)
f_V_int = [zeros(1, n_int); -gamma; zeros(1, n_int)];
% vector of volume forces
f = ASSEMBLY.vector_volume_3D(elem, coord, f_V_int, HatP, WEIGHT);

%% Input parameters for the indirect continuation
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
%% Computation of the limit load factor by the indirect continuation 

fprintf('\n Indirect continuation method for the LL method\n');
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
title('Indirect continuation method for the LL method', 'Interpreter', 'latex')
xlabel('control variable - $\omega$', 'Interpreter', 'latex');
ylabel('load factor - $t$', 'Interpreter', 'latex');