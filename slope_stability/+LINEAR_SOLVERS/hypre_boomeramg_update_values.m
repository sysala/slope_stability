function setup_time = hypre_boomeramg_update_values(V, instance_id)
%HYPRE_BOOMERAMG_UPDATE_VALUES  Update matrix values and re-run AMG setup.
%
%   setup_time = LINEAR_SOLVERS.hypre_boomeramg_update_values(V, instance_id)
%
% V must have the same length and order as the original (I,J,V) triplets
% passed to hypre_boomeramg_setup_ijv.  Only the numerical values of the
% HYPRE IJ matrix are updated (the sparsity pattern is unchanged).
% BoomerAMG setup is re-executed.

if nargin < 2
    error('Usage: hypre_boomeramg_update_values(V, instance_id)');
end

if exist('hypre_boomeramg_mex', 'file') ~= 3
    error(['hypre_boomeramg_mex is not built. Run ', ...
        'LINEAR_SOLVERS.build_hypre_boomeramg_mex from slope_stability first.']);
end

setup_time = hypre_boomeramg_mex('update_values', V(:), instance_id);
end
