function [alpha, timing] = damping_ALG5(it_damp_max, U_it, lambda_it, d_U, d_l, f, criterion, Q, constitutive_matrix_builder)
%--------------------------------------------------------------------------
% damping_ALG5 computes the damping coefficient for ALG5.
%
% This function performs a line search to determine an appropriate damping
% factor alpha such that, when updating the displacement U and the parameter
% lambda, the new residual (measured by norm(F_alpha(Q)-f(Q))) is reduced
% relative to the current criterion.
%
% INPUTS:
%   it_damp_max             - Maximum number of damping iterations.
%   U_it                    - Current displacement field.
%   lambda_it               - Current value of lambda.
%   d_U                     - Newton increment for the displacement field.
%   d_l                     - Newton increment for lambda.
%   f                       - External load vector.
%   criterion               - Current residual criterion (norm(F(Q)-f(Q))).
%   Q                       - Logical array restricting active degrees of freedom.
%   constitutive_matrix_builder - Object that builds the constitutive force
%                                 vector F given lambda and U.
%
% OUTPUTS:
%   alpha                   - Damping factor (0 <= alpha <= 1).
%   timing                  - Struct with fields:
%                               t_build_F  - total time in build_F_all calls
%                               t_norm     - total time in norm evaluation
%                               n_iters    - number of damping iterations
%
%--------------------------------------------------------------------------
if isnan(d_l)
    alpha = 0;
    timing = struct('t_build_F', 0, 't_norm', 0, 'n_iters', 0);
    return
end

alpha = 1;
it_damp = 0;
t_build_F = 0;
t_norm = 0;

while true
    it_damp = it_damp + 1;

    % Compute updated displacement and lambda using the damping factor.
    U_alpha = U_it + alpha * d_U;
    lambda_alpha = lambda_it + alpha * d_l;

    % Compute the residual at the updated values.
    t_tmp = tic;
    F_alpha = constitutive_matrix_builder.build_F_all(lambda_alpha, U_alpha);
    t_build_F = t_build_F + toc(t_tmp);

    t_tmp = tic;
    crit_alpha = norm(F_alpha(Q) - f(Q));
    t_norm = t_norm + toc(t_tmp);

    % If the new residual is not improved (i.e., not lower than the previous
    % criterion), reduce the damping factor.
    if crit_alpha >= criterion
        alpha = alpha / 2;
    else
        break  % Accept the damping factor.
    end

    % Terminate if maximum damping iterations are reached.
    if it_damp >= it_damp_max
        alpha = 0;
        break
    end
end % while

timing.t_build_F = t_build_F;
timing.t_norm    = t_norm;
timing.n_iters   = it_damp;

end % function
