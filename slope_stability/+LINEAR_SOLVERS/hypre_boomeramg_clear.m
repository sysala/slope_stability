function hypre_boomeramg_clear(instance_id)
%HYPRE_BOOMERAMG_CLEAR Clear one or all persistent HYPRE instances.
%
%   LINEAR_SOLVERS.hypre_boomeramg_clear(instance_id)
%   LINEAR_SOLVERS.hypre_boomeramg_clear()              % clear all

if ~LINEAR_SOLVERS.is_hypre_mex_available()
    return;
end

if nargin < 1
    LINEAR_SOLVERS.hypre_boomeramg_mex('clear');
else
    LINEAR_SOLVERS.hypre_boomeramg_mex('clear', instance_id);
end

end
