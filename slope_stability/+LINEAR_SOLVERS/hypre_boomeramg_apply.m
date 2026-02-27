function y = hypre_boomeramg_apply(r, instance_id)
%HYPRE_BOOMERAMG_APPLY Apply persistent BoomerAMG preconditioner instance.
%
%   y = LINEAR_SOLVERS.hypre_boomeramg_apply(r, instance_id)
%
% Inputs:
%   r           - right-hand side vector (length n).
%   instance_id - id used in LINEAR_SOLVERS.hypre_boomeramg_setup.
%
% Output:
%   y - result of one preconditioner application.

if nargin < 2
    error('Usage: hypre_boomeramg_apply(r, instance_id)');
end
if exist('hypre_boomeramg_mex', 'file') ~= 3
    error(['hypre_boomeramg_mex is not built. Run ', ...
           'LINEAR_SOLVERS.build_hypre_boomeramg_mex from slope_stability first.']);
end

y = hypre_boomeramg_mex('apply', r, instance_id);

end
