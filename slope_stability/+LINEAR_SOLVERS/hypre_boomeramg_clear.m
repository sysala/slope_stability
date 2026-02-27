function hypre_boomeramg_clear(instance_id)
%HYPRE_BOOMERAMG_CLEAR Clear one or all persistent HYPRE instances.
%
%   LINEAR_SOLVERS.hypre_boomeramg_clear(instance_id)
%   LINEAR_SOLVERS.hypre_boomeramg_clear()              % clear all

if exist('hypre_boomeramg_mex', 'file') ~= 3
    return;
end

if nargin < 1
    hypre_boomeramg_mex('clear');
else
    hypre_boomeramg_mex('clear', instance_id);
end

end
