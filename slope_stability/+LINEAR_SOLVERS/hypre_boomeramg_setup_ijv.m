function info = hypre_boomeramg_setup_ijv(I, J, V, n, null_space, opts, instance_id)
%HYPRE_BOOMERAMG_SETUP_IJV  Setup BoomerAMG from COO triplets (I,J,V).
%
%   info = LINEAR_SOLVERS.hypre_boomeramg_setup_ijv(I, J, V, n, null_space, opts, instance_id)
%
% Like hypre_boomeramg_setup but accepts COO triplets instead of a sparse
% matrix.  I, J are 1-based integer-valued double vectors; V is the values
% vector.  n is the matrix dimension.
%
% The stored (I,J) pattern is remembered so that subsequent calls to
% hypre_boomeramg_update_values can accept only the new V vector.

if nargin < 7
    error('Usage: hypre_boomeramg_setup_ijv(I, J, V, n, null_space, opts, instance_id)');
end
if isempty(opts), opts = struct(); end

if ~LINEAR_SOLVERS.is_hypre_mex_available()
    error(['hypre_boomeramg_mex is not built. Run ', ...
        'LINEAR_SOLVERS.build_hypre_boomeramg_mex from slope_stability first.']);
end

info = LINEAR_SOLVERS.hypre_boomeramg_mex('setup_ijv', I(:), J(:), V(:), double(n), ...
    null_space, opts, instance_id);
end
