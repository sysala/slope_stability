% Exported from slope_stability_3D_hetero_seepage_SSR_comsol_demo.ipynb
% for automated testing of sparsersb + HYPRE IJV path.

close all; clearvars; clc;

% Load sparsersb for multithreaded sparse matrix-vector products
pkg load sparsersb;

fprintf('Working directory: %s\n', pwd);

%% 1) Main Input Data
elem_type='P2';
Davis_type='B';

mat_props = [15, 30,  0, 10000, 0.33, 19, 19;
    15, 38,  0, 50000, 0.30, 22, 22;
    10, 35,  0, 50000, 0.30, 21, 21;
    18, 32,  0, 20000, 0.33, 20, 20;
    ];

k = 1.0;

x1 = 15;  x2 = 10;  x3 = 15;
y1 = 10;  y2 = 10;  z  = 5;
N_h = 1;

%% 2) Reference Element Data and Mesh Load
[Xi, WF] = ASSEMBLY.quadrature_volume_3D(elem_type);
[HatP, DHatP1, DHatP2, DHatP3] = ASSEMBLY.local_basis_volume_3D(elem_type, Xi);

file_path = 'meshes/comsol_mesh.h5';
[coord, elem, surf, Q, material, triangle_labels] = MESH.load_mesh_P2(file_path, 1);

n_n = size(coord,2);
n_unknown = length(coord(Q));
n_e = size(elem,2);
n_q = length(WF);
n_int = n_e * n_q;

fprintf('\nMesh data:');
fprintf('  nodes=%d  unknowns=%d  elements=%d  int_pts=%d\n', n_n, n_unknown, n_e, n_int);

material_identifier = zeros(1, n_e);

%% 3) Seepage Problem
conduct0 = SEEPAGE.heter_conduct(material_identifier, n_q, k);
grho = 9.81;
[Q_w, pw_D] = MESH.seepage_boundary_3D_hetero_comsol(coord, surf, triangle_labels, grho);

[pw, grad_p, mater_sat] = SEEPAGE.seepage_problem_3D( ...
    coord, elem, Q_w, pw_D, grho, conduct0, HatP, DHatP1, DHatP2, DHatP3, WF);

mater_sat_ext = repmat(mater_sat, n_q, 1);
saturation = mater_sat_ext(:);

%% 4) Mechanical Material Fields and Assembly
fields = {'c0', 'phi', 'psi', 'young', 'poisson', 'gamma_sat', 'gamma_unsat'};
materials = cellfun(@(x) cell2struct(num2cell(x), fields, 2), num2cell(mat_props, 2), 'UniformOutput', false);

[c0, phi, psi, shear, bulk, lame, gamma] = ...
    ASSEMBLY.heterogenous_materials(material, saturation, n_q, materials);

[K_elast, B, WEIGHT] = ASSEMBLY.elastic_stiffness_matrix_3D( ...
    elem, coord, shear, bulk, DHatP1, DHatP2, DHatP3, WF);

f_V_int = [zeros(1, n_int); -gamma; zeros(1, n_int)];
f_V = ASSEMBLY.vector_volume_3D(elem, coord, f_V_int, HatP, WEIGHT);

%% 5) Continuation, Newton, and Linear Solver Parameters
lambda_init = 1.0;
d_lambda_init = 0.1;
d_lambda_min = 1e-5;
d_lambda_diff_scaled_min = 0.005;
omega_max_stop = 340700000;
step_max = 2;

it_newt_max = 50;
it_damp_max = 10;
tol = 1e-4;
r_min = 1e-4;

agmg_folder = "agmg";
solver_type = 'DFGMRES_HYPRE_BOOMERAMG';

linear_solver_tolerance = 1e-1;
linear_solver_maxit = 100;
deflation_basis_tolerance = 1e-3;
linear_solver_printing = 0;

boomeramg_opts = struct('threads', 16, 'print_level', 0, ...
    'use_as_preconditioner', true);

linear_system_solver = LINEAR_SOLVERS.set_linear_solver(agmg_folder, solver_type, ...
    linear_solver_tolerance, linear_solver_maxit, deflation_basis_tolerance, ...
    linear_solver_printing, Q, coord, boomeramg_opts);

n_strain = 6;
constitutive_matrix_builder = CONSTITUTIVE_PROBLEM.CONSTITUTIVE( ...
    B, c0, phi, psi, Davis_type, shear, bulk, lame, WEIGHT, n_strain, n_int, 3);

%% 6) Run SSR Continuation
fprintf('\n Indirect continuation method\n');
tic;
[U3, lambda_hist3, omega_hist3, Umax_hist3, stats] = CONTINUATION.SSR_indirect_continuation( ...
    lambda_init, d_lambda_init, d_lambda_min, d_lambda_diff_scaled_min, step_max, ...
    omega_max_stop, it_newt_max, it_damp_max, tol, r_min, K_elast, Q, f_V, ...
    constitutive_matrix_builder, linear_system_solver.copy());
time_run = toc;
fprintf('Running_time = %f \n', time_run);

if ~isempty(strfind(upper(char(solver_type)), 'BOOMERAMG'))
    LINEAR_SOLVERS.hypre_boomeramg_clear();
end

%% Newton Profiler Summary
if exist('stats', 'var') && isfield(stats, 'newton_timing')
    nt = stats.newton_timing;
    labels = fieldnames(nt);
    times  = cellfun(@(f) nt.(f), labels);
    total  = sum(times);

    [times, idx] = sort(times, 'descend');
    labels = labels(idx);
    pcts   = 100 * times / total;

    desc_map = struct( ...
        'build_F_and_DS',      'build_F_and_DS_all     (constitutive F + DS)', ...
        'solve_V',             'linear_solver.solve     (V solve)',           ...
        'A_orthogonalize',     'linear_solver.A_orthogonalize',              ...
        'damping',             'NEWTON.damping_ALG5     (line search)',       ...
        'setup_preconditioner','setup/update preconditioner (HYPRE IJV)',     ...
        'build_F_eps',         'build_F_all             (F_eps for diff)',    ...
        'solve_W',             'linear_solver.solve     (W solve)',           ...
        'K_r_assembly',        'K_r sparsersb assembly  (BtDB prebuilt)',    ...
        'expand_deflation_W',  'expand_deflation_basis  (W)',                ...
        'expand_deflation_V',  'expand_deflation_basis  (V)');

    fprintf('\n============================================================\n');
    fprintf('  Newton Profiler  (newton_ind_SSR accumulated)\n');
    fprintf('============================================================\n');
    fprintf('  %-6s  %-8s  %s\n', 'Time', '  %', 'Operation');
    fprintf('  %-6s  %-8s  %s\n', '------', '------', '------------------------------');
    for k = 1:numel(labels)
        lbl = labels{k};
        if isfield(desc_map, lbl)
            desc = desc_map.(lbl);
        else
            desc = lbl;
        end
        fprintf('  %6.1fs  %5.1f%%   %s\n', times(k), pcts(k), desc);
    end
    fprintf('  %-6s  %-8s  %s\n', '------', '------', '------------------------------');
    fprintf('  %6.1fs  %5.1f%%   %s\n', total, 100.0, 'TOTAL');
    fprintf('============================================================\n');
end

fprintf('\nTest completed successfully.\n');
