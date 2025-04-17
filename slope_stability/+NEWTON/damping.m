function alpha = damping(it_damp_max, U_it, dU, F, f, constitutive_matrix_builder)
%--------------------------------------------------------------------------
% damping computes the damping coefficient for Newton's iteration.
%
% This function performs a line search to find an appropriate damping factor
% alpha that yields a descent direction. It uses a bisection strategy on the 
% interval [alpha_min, alpha_max] to ensure that the directional derivative 
% of the residual (F(U)-f) along the increment dU is negative.
%
% INPUTS:
%   it_damp_max             - Maximum number of damping iterations.
%   U_it                    - Current displacement field.
%   dU                      - Newton increment for the displacement.
%   F                       - Current internal force vector F(U_it).
%   f                       - External load vector.
%   constitutive_matrix_builder - Object to build the constitutive force vector,
%                                 used here to recompute F at updated U.
%
% OUTPUT:
%   alpha                   - Damping factor (0 <= alpha <= 1). A value of 0
%                             indicates failure to obtain a descent direction.
%
%--------------------------------------------------------------------------
% Check initial descent condition.
decrease = (F(:) - f(:))' * dU(:);
if isnan(decrease) || (norm(dU(:)) == Inf) || (decrease >= 0)
    alpha = 0;
    return
end

alpha = 1; 
alpha_min = 0; 
alpha_max = 1;
it_damp = 0;

while it_damp < it_damp_max
    it_damp = it_damp + 1;
    % Compute the updated displacement.
    U_alpha = U_it + alpha * dU;
    % Recompute the internal force at the updated displacement.
    F_alpha = constitutive_matrix_builder.build_F_reduced(U_alpha);
    % Compute the directional derivative at the new iterate.
    decrease = (F_alpha(:) - f(:))' * dU(:);
    if decrease < 0
        if alpha == 1
            break;  % Accept full step if descent is observed.
        end
        alpha_min = alpha;  % Increase lower bound.
    else
        alpha_max = alpha;  % Decrease upper bound.
    end
    % Bisection update of the damping factor.
    alpha = (alpha_min + alpha_max) / 2;
end

end % function
