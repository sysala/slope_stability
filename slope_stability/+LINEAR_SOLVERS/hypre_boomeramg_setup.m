function info = hypre_boomeramg_setup(A, null_space, opts, instance_id)
%HYPRE_BOOMERAMG_SETUP Setup persistent BoomerAMG preconditioner instance.
%
%   info = LINEAR_SOLVERS.hypre_boomeramg_setup(A, null_space, opts, instance_id)
%
% Inputs:
%   A           - real sparse square matrix (double).
%   null_space  - dense matrix [n x k] with near-null-space vectors
%                 (e.g., rigid-body modes), or [].
%   opts        - struct with BoomerAMG options, or [] for defaults.
%   instance_id - numeric scalar or char id. If the id already exists,
%                 it is destroyed and reinitialized with current inputs.
%
% Persistent behavior:
%   Instances remain alive while MATLAB runs (until explicitly cleared or
%   mex is unloaded). Multiple instances are supported via instance_id.
%
% Main opts fields supported (all optional):
%   threads, print_level, max_levels,
%   coarsen_type, agg_num_levels, interp_type, p_max_elmts,
%   strong_threshold, strong_threshold_r,
%   relax_type, relax_sweeps, relax_coarse_type, cycle_type,
%   num_functions, dof_func,
%   nodal, nodal_diag, nodal_levels, keep_same_sign,
%   interp_vec_variant, interp_vec_qmax, interp_vec_abs_qtrunc,
%   smooth_interp_vectors, interp_refine,
%   trunc_factor, agg_interp_type, agg_trunc_factor, agg_p_max_elmts,
%   max_coarse_size, min_coarse_size,
%   use_as_preconditioner, max_iter, tol.
%
% Default configuration when opts = [] (elasticity tuned):
%   threads = 16
%   use_as_preconditioner = true  -> one BoomerAMG apply (max_iter=1, tol=0)
%   coarsen_type = 8, agg_num_levels = 0, interp_type = 17
%   strong_threshold = 0.5, p_max_elmts = 4
%   relax_type = 8, relax_sweeps = 1, relax_coarse_type = 8, cycle_type = 1
%   num_functions = 3 if null_space provided and n mod 3 == 0, else 1
%   nodal = 4, nodal_diag = 1
%   interp_vec_variant = 2, interp_vec_qmax = 4, smooth_interp_vectors = 1
%
% Output:
%   info struct fields:
%     instance_id, n, nnz, num_functions, num_interp_vectors,
%     setup_time_seconds
%
% Companion calls:
%   y = LINEAR_SOLVERS.hypre_boomeramg_apply(r, instance_id)
%   LINEAR_SOLVERS.hypre_boomeramg_clear(instance_id)   % or clear all
%
% Notes:
%   This API is designed for preconditioner application (one cycle).
%   For custom solver behavior set use_as_preconditioner=false and pass
%   max_iter/tol explicitly.

if nargin < 4
    error('Usage: hypre_boomeramg_setup(A, null_space, opts, instance_id)');
end
if isempty(opts)
    opts = struct();
end
if exist('hypre_boomeramg_mex', 'file') ~= 3
    error(['hypre_boomeramg_mex is not built. Run ', ...
           'LINEAR_SOLVERS.build_hypre_boomeramg_mex from slope_stability first.']);
end

info = hypre_boomeramg_mex('setup', A, null_space, opts, instance_id);

end
