function [shur_prec_handle] = diag_prec_HYPRE_BOOMERAMG(A, null_space, opts, instance_id)
% diag_prec_HYPRE_BOOMERAMG builds a BoomerAMG preconditioner handle via MEX.
%
% Inputs:
%   A           - Sparse system matrix.
%   null_space  - Near-null-space vectors (n x k), or [].
%   opts        - HYPRE options struct passed to hypre_boomeramg_setup.
%   instance_id - Persistent instance id for the HYPRE MEX backend.
%
% Output:
%   shur_prec_handle - Function handle y = M^{-1}x for use in DFGMRES.

if nargin < 2 || isempty(null_space)
    null_space = [];
end
if nargin < 3 || isempty(opts)
    opts = struct();
end
if nargin < 4 || isempty(instance_id)
    instance_id = 'hypre_boomeramg_default';
end

LINEAR_SOLVERS.hypre_boomeramg_setup(A, null_space, opts, instance_id);
shur_prec_handle = @(x) LINEAR_SOLVERS.hypre_boomeramg_apply(x, instance_id);

end
