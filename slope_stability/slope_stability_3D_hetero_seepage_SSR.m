%%  Stability of a homogeneous slope with seepage in 3D
% =========================================================================
%
%  This program solves a 3D slope stability problem by the modified shear
%  strength reduction (SSR) method described in (Sysala et al., CAS 2025). 
%  The Mohr-Coulomb yield criterion, 3 Davis approaches (denoted by A, B, C),
%  standard finite elements (P1 or P2) and meshes with different densities 
%  are considered. For P2 elements, the 11-point Gauss quadrature is used.
%  To find the safety factor of the SSR method, two continuation techniques
%  are available: direct and indirect. A bechmark problem on a homogeneous 
%  slope with unconfined seepage is considered. It is possible to change
%  geometrical parameters and mesh density.
%
% ======================================================================
%

%% The main input data

% elem_type - type of finite elements; available choices: 'P1', 'P2'
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
mat_props = [6, 45, 0, 40000, 0.3, 20, 20];  

% Hydraulic conductivity for each subdomain [m/s]
k = 1.0 ;

% Geometrical parameters (choose only integers except beta)
  x1 = 15 ;        % length of the body in front of the slope
  x2 = 10 ;        % length of the the slope in x-direction
  x3 = 15 ;        % length of the body behind the slope
  y1 = 10 ;        % hight of the body below the slope
  y2 = 10 ;        % height of the slope
  z = 5;           % length of the body in z-direction
  
% Mesh density parameter (choose only integer!)
  N_h  = 1;        % h=1/N_h - discretization parameter

  
%% Data from the reference element
% quadrature points and weights for volume integration
[Xi, WF] = ASSEMBLY.quadrature_volume_3D(elem_type);
% local basis functions and their derivatives
[HatP, DHatP1, DHatP2, DHatP3] = ASSEMBLY.local_basis_volume_3D(elem_type, Xi);

%% Creation of the uniform finite element mesh
% switch(elem_type)
%     case 'P1'
%         [coord, elem, surf, Q] = MESH.mesh_P1_3D(N_h, x1, x2, x3, y1, y2, z);
%         fprintf('P1 elements: \n')
%     case 'P2'
%         [coord, elem, surf, Q] = MESH.mesh_P2_3D(N_h, x1, x2, x3, y1, y2, z);
%         fprintf('P2 elements: \n')
%     otherwise
%         error('Bad choice of element type');
% end
[coord, elem, surf, Q, material, triangle_labels] = MESH.load_mesh_gmsh_waterlevels("meshes/slope_with_waterlevels.h5");


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

% The array material_identifier for a homogeneous body
material_identifier = zeros(1,n_e);

%% Computation of porous water pressure

% Hydraulic conductivity ateach integration point
conduct0=SEEPAGE.heter_conduct(material_identifier,n_q,k); 

% specific weight of water in kPa
grho=9.81;

% Dirichlet boundary conditions for pressure (problem dependent)
% Q_w=true(1,n_n);
% Q_w(coord(1,:)<=0.001)=0;
% Q_w(coord(1,:)>=x1+x2+x3-0.001)=0;
% Q_w(coord(2,:)>=y1+y2-0.001)=0;
% Q_w((coord(2,:)>=y1-0.001)&(coord(1,:)>=x1+x2-0.001))=0;
% Q_w((coord(2,:)>=y1-0.001)&(coord(2,:)>=-(y2/x2)*coord(1,:)+y1+y2*(1+x1/x2)-0.001))=0;  

% Nonhomogeneous part of the pressure (problem dependent)
% y21=2;            % height of the water level next to the slope
% y22=6;            % difference between water levels on oposite slope sides
% y23=2;            % height of the slope above underground water level
% pw_D=zeros(1,n_n);
% x_bar=x1+(1-y21/y2)*x2;
% part1=(coord(1,:)<x_bar)&(coord(2,:)<=-(y22/x_bar)*coord(1,:)+y1+y21+y22);
% part2=coord(1,:)>=x_bar;
% pw_D(part1)=grho*((y22/x_bar)*(x_bar-coord(1,part1))+y1+y21-coord(2,part1));
% pw_D(part2)=grho*(y1+y21-coord(2,part2)); 

[Q_w, pw_D] = MESH.darcy_boundary_3D_hetero(coord, surf, triangle_labels, grho);

% Computation on the pore pressure and its gradient
[pw, grad_p, mater_sat]=SEEPAGE.seepage_problem_3D...
          (coord,elem,Q_w,pw_D,grho,conduct0,HatP,DHatP1,DHatP2,DHatP3,WF);

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
      ASSEMBLY.heterogenous_materials(material*0, saturation, n_q, materials);

%% Assembling for mechanics

% Assembling of the elastic stiffness matrix
[K_elast, B, WEIGHT] = ASSEMBLY.elastic_stiffness_matrix_3D(elem, coord, shear, bulk, DHatP1, DHatP2, DHatP3, WF);

% Assemble the vector of volume forces.
% Volume forces at integration points, size (f_V_int) = (3, n_int)
f_V_int = [-grad_p(1,:); -grad_p(2,:)-gamma; -grad_p(3,:)];
% f_V_int = [zeros(1, n_int); -gamma; zeros(1, n_int)];
% Compute the vector of volume forces.
f_V = ASSEMBLY.vector_volume_3D(elem, coord, f_V_int, HatP, WEIGHT);


%% Input parameters for the continuation methods

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
%% Computation of the factor of safety for the SSR method

direct_on = 0; % Use direct continuation method.
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
    % VIZ.draw_saturation_2D(coord,elem,mater_sat);
    VIZ.draw_mesh_3D(coord,surf);
    VIZ.draw_quantity_3D_old(coord,surf,zeros(size(coord)),pw,x1,x2,x3,y1,y2,z);
    % VIZ.plot_pore_pressure_3D(pw,coord,elem);
    VIZ.plot_displacements_3D(U3, coord, elem);
    VIZ.plot_deviatoric_strain_3D(U3, coord, elem, B);
    % Visualization of the curve: omega -> lambda for direct continuation.
    figure; hold on; box on; grid on;
    plot(omega_hist2, lambda_hist2, '-o');
    title('Direct continuation method', 'Interpreter', 'latex')
    xlabel('variable - $\xi$', 'Interpreter', 'latex');
    ylabel('strength reduction factor - $\lambda$', 'Interpreter', 'latex');
end

%% Postprocessing - visualization of selected results for indirect continuation
if indirect_on
    % VIZ.draw_saturation_2D(coord,elem,mater_sat);
    VIZ.draw_mesh_3D(coord,surf);
    VIZ.draw_quantity_3D_old(coord,surf,zeros(size(coord)),pw,x1,x2,x3,y1,y2,z);
    % VIZ.plot_pore_pressure_3D(pw,coord,elem);
    VIZ.plot_displacements_3D(U3, coord, elem);
    VIZ.plot_deviatoric_strain_3D(U3, coord, elem, B);
    % Visualization of the curve: omega -> lambda for indirect continuation.
    figure; hold on; box on; grid on;
    plot(omega_hist3, lambda_hist3, '-o');
    title('Indirect continuation method', 'Interpreter', 'latex')
    xlabel('control variable - $\omega$', 'Interpreter', 'latex');
    ylabel('strength reduction factor - $\lambda$', 'Interpreter', 'latex');
end
