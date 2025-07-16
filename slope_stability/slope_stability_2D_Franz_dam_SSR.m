%%  Stability of a high heterogeneous embankment dam via the SSR methods
% =========================================================================
%
%  This program solves a 2D slope stability problem by the modified shear
%  strength reduction (SSR) method described in (Sysala et al., CAS 2025). 
%  The Mohr-Coulomb yield criterion, 3 Davis approaches (denoted by A, B, C),
%  standard finite elements (P1, P2 or P4 elements) and meshes
%  with different densities are considered. For P2 elements, the 7-point 
%  Gauss quadrature is used. To find the safety factor of the SSR method, 
%  two continuation techniques are available: direct and indirect. A
%  bechmark problem on a high heterogeneous embankment dam with unconfined 
%  seepage is considered, see (Sysala et al., CAS 2023).
%
% ======================================================================
%

%% The main input data

% elem_type - type of finite elements; available choices: 'P1', 'P2', 'P4'
elem_type='P2';

% Davis_type - choice of Davis' approach; available choices: 'A','B','C'
Davis_type='B';

% Mechanical parameters for each subdomain. In the following table, we
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
mat_props = ...
   [50.0, 42.00, 42.00, 16000, 0.4, 27.0, 27.0;  % Subdomain 1 - Zone X
     0.5, 41.00,  0.00, 16000, 0.4, 22.0, 22.0;  % Subdomain 2 - Zone RF
     0.5, 40.00,  0.00, 16000, 0.4, 23.0, 23.0;  % Subdomain 3 - Zone C2
     0.0, 38.00,  0.00, 16000, 0.4, 23.0, 23.0;  % Subdomain 4 - Zone B
    10.0, 25.00,  0.00, 16000, 0.4, 22.0, 22.0;  % Subdomain 5 - Zone A
     0.5, 41.00,  0.00, 16000, 0.4, 22.0, 22.0;  % Subdomain 6 - Zone RF
    75.0, 42.00, 42.00, 16000, 0.4, 27.0, 27.0;  % Subdomain 7 - Zone E
     0.0, 38.00,  0.00, 16000, 0.4, 23.0, 23.0;  % Subdomain 8 - Zone B
     0.5, 41.00,  0.00, 16000, 0.4, 23.0, 23.0;  % Subdomain 9 - Zone C1
     0.5, 41.00,  0.00, 16000, 0.4, 21.0, 21.0]; % Subdomain 10 - Zone D  

%  Hydraulic conductivity for each subdomain [m/s]
  k = [1.0e-5   % Subdomain 1 - Zone X
       1.0e-4   % Subdomain 2 - Zone RF
       5.0e-5   % Subdomain 3 - Zone C2
       5.0e-6   % Subdomain 4 - Zone B
       1.0e-9   % Subdomain 5 - Zone A
       1.0e-4   % Subdomain 6 - Zone RF
       2.0e-9   % Subdomain 7 - Zone E
       5.0e-6   % Subdomain 8 - Zone B
       5.0e-5   % Subdomain 9 - Zone C1     
       5.0e-4]; % Subdomain 10 - Zone D            

%% Data from the reference element
% quadrature points and weights for volume integration
[Xi, WF] = ASSEMBLY.quadrature_volume_2D(elem_type);
% local basis functions and their derivatives
[HatP,DHatP1,DHatP2] = ASSEMBLY.local_basis_volume_2D(elem_type, Xi);

%% Creation/loading of the finite element mesh
[coord, elem, Q, material_identifier, surf] = MESH.load_mesh_Franz_dam(elem_type, 'meshes/Franz_dam/');
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

%% Computation of porous water pressure

% Hydraulic conductivity ateach integration point
conduct0=SEEPAGE.heter_conduct(material_identifier,n_q,k); 

% specific weight of water in kPa
grho=9.81;

% Dirichlet boundary conditions for pressure (problem dependent)
Q_D=(coord(2,surf(1,:))+coord(2,surf(2,:)))/2>=-400+1e-1;
Q_w=true(1,n_n);
Q_w(unique(surf(:,Q_D)))=0;

% Nonhomogeneous part of the pressure (problem dependent)
x1=-82.5; y1=-50; 
x2=172.5; y2=-112; 
pw_D=zeros(1,n_n);
part1=(coord(1,:)<x1+1e-9)&(coord(2,:)<y1);
part2=(coord(1,:)>=x1+1e-9)&(coord(1,:)<x2+1e-9)&(coord(2,:)<((y2-y1)/(x2-x1))*(coord(1,:)-x1)+y1);
part3=(coord(1,:)>=x2+1e-9)&(coord(2,:)<y2);
pw_D(part1)=grho*(y1-coord(2,part1));
pw_D(part2)=grho*(((y2-y1)/(x2-x1))*(coord(1,part2)-x1)+y1-coord(2,part2));
pw_D(part3)=grho*(y2-coord(2,part3));  

% Computation on the pore pressure and its gradient
[pw, grad_p, mater_sat]=SEEPAGE.seepage_problem_2D...
                         (coord,elem,Q_w,pw_D,grho,conduct0,HatP,DHatP1,DHatP2,WF);

% Saturation - a prescribed logical array indicating integration points 
%              where the body is saturated. If gamma_sat and gamma_unsat 
%              are the same, set saturation=true(1,n_int). Otherwise,
%              this logical array is derived from the phreatic surface.
mater_sat_ext=repmat(mater_sat,n_q,1);
saturation=mater_sat_ext(:);

%% Mechanical material Parameters at Integration Points

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

% Material parameters at integration points.
[c0, phi, psi, shear, bulk, lame, gamma] = ...
      ASSEMBLY.heterogenous_materials(material_identifier, saturation, n_q, materials);

%% Assembling for mechanics

% Assembling of the elastic stiffness matrix
[K_elast,B,WEIGHT]=ASSEMBLY.elastic_stiffness_matrix_2D(elem,coord,...
    DHatP1,DHatP2,WF,shear,lame);

% volume forces at integration points, size(f_V_int)=(2,n_int)
% grad_p=zeros(2,n_int);
f_V_int = [-grad_p(1,:);-grad_p(2,:)-gamma] ;
% vector of volume forces
f_V=ASSEMBLY.vector_volume_2D(elem,coord,f_V_int,HatP,WEIGHT);


%% Input parameters for the continuation methods

lambda_init = 0.7;              % Initial lower bound of lambda
d_lambda_init = 0.1;            % Initial increment of lambda
d_lambda_min = 1e-5;            % Minimal increment of lambda
d_lambda_diff_scaled_min = 0.001;% Minimal rate of increment of lambda
omega_max_stop = 1e10;           % Maximum omega, then stop
step_max = 100;                 % Maximum number of continuation steps

%% Input parameters for Newton's solvers
it_newt_max = 50;               % Number of Newton's iterations
it_damp_max = 10;               % Number of iterations within line search
tol = 1e-5;                     % Relative tolerance for Newton's solvers
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

direct_on = 1; % Use direct continuation method.
indirect_on = 1; % Use indirect continuation method.

if direct_on  % Direct continuation method.
    fprintf('\n Direct continuation method\n');
    tic;
    [U2, lambda_hist2, omega_hist2, Umax_hist2] = CONTINUATION.SSR_direct_continuation(...
        lambda_init, d_lambda_init, d_lambda_min, d_lambda_diff_scaled_min, step_max, ...
        it_newt_max, it_damp_max, tol, r_min, K_elast, Q, f_V, ...
        constitutive_matrix_builder, linear_system_solver.copy());
    time_run = toc;
    fprintf("Running_time = %f \n", time_run);
end
if indirect_on     % Indirect continuation method.
    fprintf('\n Indirect continuation method\n');
    tic;
    [U3, lambda_hist3, omega_hist3, Umax_hist3] = CONTINUATION.SSR_indirect_continuation(...
        lambda_init, d_lambda_init, d_lambda_min, d_lambda_diff_scaled_min, step_max, ...
        omega_max_stop, it_newt_max, it_damp_max, tol, r_min, K_elast, Q, f_V, ...
        constitutive_matrix_builder, linear_system_solver.copy());
    time_run = toc;
    fprintf("Running_time = %f \n", time_run);
end

%% Postprocessing - visualization of selected results for direct continuation
if direct_on
    VIZ.draw_heterogeneity_Franz_dam(coord,elem,material_identifier);
    VIZ.draw_mesh_2D(coord,elem)
    VIZ.plot_pore_pressure_2D(pw,coord,elem);
    VIZ.draw_saturation_2D(coord,elem,mater_sat);
    VIZ.plot_deviatoric_strain_2D(U2,coord,elem,B);
    VIZ.plot_displacements_2D(U2,coord,elem);
    % Visualization of the curve: omega -> lambda for direct continuation.
    figure; hold on; box on; grid on;
    plot(omega_hist2, lambda_hist2, '-o');
    title('Direct continuation method', 'Interpreter', 'latex')
    xlabel('variable - $\xi$', 'Interpreter', 'latex');
    ylabel('strength reduction factor - $\lambda$', 'Interpreter', 'latex');
end

%% Postprocessing - visualization of selected results for indirect continuation
if indirect_on
    VIZ.draw_heterogeneity_Franz_dam(coord,elem,material_identifier);
    VIZ.draw_mesh_2D(coord,elem)
    VIZ.plot_pore_pressure_2D(pw,coord,elem);
    VIZ.draw_saturation_2D(coord,elem,mater_sat);
    VIZ.plot_deviatoric_strain_2D(U3,coord,elem,B);
    VIZ.plot_displacements_2D(U3,coord,elem);
    % Visualization of the curve: omega -> lambda for indirect continuation.
    figure; hold on; box on; grid on;
    plot(omega_hist3, lambda_hist3, '-o');
    title('Indirect continuation method', 'Interpreter', 'latex')
    xlabel('control variable - $\omega$', 'Interpreter', 'latex');
    ylabel('strength reduction factor - $\lambda$', 'Interpreter', 'latex');
end
