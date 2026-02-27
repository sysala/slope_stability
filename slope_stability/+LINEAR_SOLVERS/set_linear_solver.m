function [linear_system_solver] = set_linear_solver(agmg_folder, solver_type, ...
    linear_solver_tolerance, linear_solver_maxit, deflation_basis_tolerance, ...
    linear_solver_printing, Q, coord, boomeramg_opts)
%SET_LINEAR_SOLVER Factory for iterative/direct solver configuration.
%
% Existing API is preserved. Additional optional BoomerAMG inputs:
%   coord          - nodal coordinates (2/3 x n) for automatic component
%                    mapping. For 3D, elasticity near-null-space vectors
%                    are also built automatically.
%   boomeramg_opts - struct forwarded to hypre_boomeramg_setup.
%                    Optional local-only fields:
%                      instance_id       - persistent HYPRE instance id
%                      null_space        - explicit near-null-space basis
%                      center_coordinates - center coords for rotation modes

if nargin < 7
    Q = [];
end
if nargin < 8
    coord = [];
end
if nargin < 9 || isempty(boomeramg_opts)
    boomeramg_opts = struct();
end
if ~isstruct(boomeramg_opts)
    error('boomeramg_opts must be a struct or empty.');
end

solver_type = upper(string(solver_type));
is_agmg_present = LINEAR_SOLVERS.check_agmg_present(agmg_folder);

if solver_type == "DFGMRES_AGMG" && is_agmg_present == 0
    warning('AGMG not found. Switching to DFGMRES_ICHOL.');
    solver_type = "DFGMRES_ICHOL";
end

switch solver_type
    case "DIRECT"
        linear_system_solver = LINEAR_SOLVERS.DIRECT_BACKSLASH();

    case "DFGMRES_ICHOL"
        preconditioner_builder = @(A) LINEAR_SOLVERS.diag_prec_ICHOL(A, Q);
        linear_system_solver = LINEAR_SOLVERS.DFGMRES(preconditioner_builder, ...
            linear_solver_tolerance, linear_solver_maxit, ...
            deflation_basis_tolerance, linear_solver_printing);

    case "DFGMRES_AGMG"
        preconditioner_builder = @(A) LINEAR_SOLVERS.diag_prec_AGMG(A, Q);
        linear_system_solver = LINEAR_SOLVERS.DFGMRES(preconditioner_builder, ...
            linear_solver_tolerance, linear_solver_maxit, ...
            deflation_basis_tolerance, linear_solver_printing);

    case {"DFGMRES_HYPRE_BOOMERAMG", "DFGMRES_BOOMERAMG", "BOOMERAMG", "HYPRE_BOOMERAMG"}
        preconditioner_builder = local_build_boomer_builder(Q, coord, boomeramg_opts);
        linear_system_solver = LINEAR_SOLVERS.DFGMRES(preconditioner_builder, ...
            linear_solver_tolerance, linear_solver_maxit, ...
            deflation_basis_tolerance, linear_solver_printing);

    otherwise
        error('Bad choice of solver type: %s', char(solver_type));
end

end

function preconditioner_builder = local_build_boomer_builder(Q, coord, opts_in)

if exist('hypre_boomeramg_mex', 'file') ~= 3
    this_dir = fileparts(mfilename('fullpath'));
    if isfile(fullfile(this_dir, 'build_hypre_boomeramg_mex.m'))
        LINEAR_SOLVERS.build_hypre_boomeramg_mex();
    else
        error(['hypre_boomeramg_mex is not built and ', ...
            'LINEAR_SOLVERS.build_hypre_boomeramg_mex is not available.']);
    end
end

opts = opts_in;
instance_id = local_make_instance_id();
if isfield(opts, 'instance_id')
    instance_id = char(string(opts.instance_id));
    opts = rmfield(opts, 'instance_id');
end

null_space = [];
if isfield(opts, 'null_space')
    null_space = opts.null_space;
    opts = rmfield(opts, 'null_space');
end

center_coordinates = true;
if isfield(opts, 'center_coordinates')
    center_coordinates = logical(opts.center_coordinates);
    opts = rmfield(opts, 'center_coordinates');
end

auto_dof_func = [];
if ~isempty(coord) && ~isempty(Q) && isequal(size(Q), size(coord)) && ...
        any(size(coord, 1) == [2, 3])
    n_components = size(coord, 1);

    % Build 3D rigid-body near-null-space automatically for elasticity.
    if n_components == 3 && isempty(null_space)
        null_space = LINEAR_SOLVERS.near_null_space_elasticity_3D(coord, Q, center_coordinates);
    end

    q_idx = find(Q(:));
    auto_dof_func = mod(q_idx - 1, n_components);

    if ~isfield(opts, 'num_functions')
        opts.num_functions = n_components;
    end
    if ~isfield(opts, 'dof_func')
        opts.dof_func = auto_dof_func(:);
    end
elseif ~isempty(coord) && isempty(null_space)
    warning(['BoomerAMG: coord/Q size mismatch, skipping automatic ', ...
        'component mapping and near-null-space construction.']);
end

if ~isfield(opts, 'print_level')
    opts.print_level = 0;
end
if ~isfield(opts, 'use_as_preconditioner')
    opts.use_as_preconditioner = true;
end

preconditioner_builder = @(A) local_build_boomer_prec(A, null_space, opts, auto_dof_func, instance_id);

end

function M = local_build_boomer_prec(A, null_space_template, opts_template, auto_dof_func, instance_id)

opts = opts_template;
n = size(A, 1);

null_space = [];
if ~isempty(null_space_template)
    if size(null_space_template, 1) == n
        null_space = null_space_template;
    else
        warning(['BoomerAMG: provided near-null-space size does not match matrix ', ...
            'dimension, ignoring near-null-space for this setup.']);
    end
end

if ~isfield(opts, 'dof_func') && ~isempty(auto_dof_func) && numel(auto_dof_func) == n
    opts.dof_func = auto_dof_func(:);
end

M = LINEAR_SOLVERS.diag_prec_HYPRE_BOOMERAMG(A, null_space, opts, instance_id);

end

function instance_id = local_make_instance_id()
[~, token] = fileparts(tempname);
instance_id = ['hypre_boomeramg_' token];
end
