function [U, lambda_hist, omega_hist, Umax_hist, work_hist] = SSR_direct_continuation(...
    lambda_init, d_lambda_init, d_lambda_min, d_lambda_diff_scaled_min, step_max, ...
    it_newt_max, it_damp_max, tol, r_min, K_elast, Q, f, ...
    constitutive_matrix_builder, linear_system_solver)
%--------------------------------------------------------------------------
% SSR_direct_continuation performs a direct continuation method to estimate 
% the safety factor lambda* and to construct the curve between the strength 
% reduction parameter lambda and the parameter omega = -G'(lambda).
%
% This method incrementally updates lambda (the strength reduction parameter)
% and computes the corresponding displacement field U and omega (work measure)
% by solving the nonlinear system:
%
%         F_lambda(U) = f,
%
% with the additional constraint that omega = -G'(lambda).
%
% The algorithm uses an initial phase to determine starting values and then 
% proceeds with continuation steps. At each step, if the solver fails or if 
% the computed omega does not increase, the increment in lambda is reduced.
%
% Alternative methods (e.g., bisection or tangent methods w.r.t. lambda) are 
% provided as commented-out code.
%
% INPUTS:
%   lambda_init              - Initial value of the parameter lambda.
%   d_lambda_init            - Initial increment of lambda.
%   d_lambda_min             - Minimal allowed increment of lambda.
%   d_lambda_diff_scaled_min - Minimal scaled increase of lambda (stopping criterion).
%   step_max                 - Maximum number of continuation steps.
%   it_newt_max              - Maximum number of Newton iterations per step.
%   it_damp_max              - Maximum number of damping iterations in Newton's method.
%   tol                      - Relative tolerance for Newton's solver.
%   r_min                    - Basic regularization parameter for the stiffness matrix.
%   K_elast                  - Elastic stiffness matrix.
%   Q                        - Logical array restricting nodes with Dirichlet boundary cond.
%   f                        - External load vector.
%   constitutive_matrix_builder - Object that builds the constitutive force vector 
%                                 and tangent operator.
%   linear_system_solver     - Linear solver object (with deflation capability).
%
% OUTPUTS:
%   U           - Displacement field corresponding to the computed lambda*.
%   lambda_hist - Array of lambda values over continuation steps.
%   omega_hist  - Array of omega values over continuation steps.
%   Umax_hist   - Array of the maximum nodal displacement per continuation step.
%
%--------------------------------------------------------------------------
%
% Initialization
%
lambda_hist = zeros(1, 1000);    % History of lambda values.
omega_hist = zeros(1, 1000);     % History of omega values.
eps = tol * 100;                  % Precision for numerical derivatives.
Umax_hist = zeros(1, 1000);      % History of maximum displacement values.
work_hist = zeros(1, 1000);

% First two steps of the continuation method.
[U, omega_old, omega, lambda, U_old] = CONTINUATION.init_phase_SSR_direct_continuation(...
    lambda_init, d_lambda_init, d_lambda_min, ...
    it_newt_max, it_damp_max, tol, eps, r_min, K_elast, Q, f, ...
    constitutive_matrix_builder, linear_system_solver);
linear_system_solver.expand_deflation_basis(U(Q));

% Storage of the computed values.
omega_hist(1) = omega_old;
lambda_hist(1) = lambda_init;
omega_hist(2) = omega;
lambda_hist(2) = lambda;
Umax_hist(1) = max(sqrt(sum(U.^2, 1)));
Umax_hist(2) = max(sqrt(sum(U.^2, 1)));
work_hist(1) = dot(U_old(:), f(:));
work_hist(2) = dot(U(:), f(:));

% Other initial data.
d_omega = omega - omega_old;     % Current increment of omega.
d_lambda = lambda - lambda_init;   % Current increment of lambda.

%
% Continuation loop over the parameter lambda.
%
step = 2;     % Current number of continuation steps.
while true
    fprintf('\n');
    fprintf(' Step = %d  ', step + 1);
    fprintf('\n');
    
    % Update of the parameter lambda.
    lambda_it = lambda + d_lambda;
    
    % U_ini = d_omega * (U - U_old) / (omega - omega_old) + U;
    % Computation of U and omega for given lambda_it.
    [U_it, omega_it, flag] = CONTINUATION.omega_SSR_direct_continuation(...
        lambda_it, U, eps, d_lambda, ...
        it_newt_max, it_damp_max, tol, r_min, K_elast, Q, f, ...
        constitutive_matrix_builder, linear_system_solver.copy());
    
    % Test on convergence:
    d_omega_test = omega_it - omega_hist(step);
    if (flag == 1) || (d_omega_test < 0)
        % The solver was not successful - decrease d_lambda.
        d_lambda = d_lambda / 2;
    % elseif d_omega_test > 2*d_omega
    %     % too large increment of omega
    %     %   U = (U + U_it) / 2;
    %     d_lambda = d_lambda / 2;
    else  % The solver was successful.
        U = U_it;
        omega = omega_it;
        lambda = lambda_it;
        step = step + 1;
        lambda_hist(step) = lambda;
        omega_hist(step) = omega;
        Umax_hist(step) = max(sqrt(sum(U.^2, 1)));
        work_hist(step) = dot(U(:), f(:));
        linear_system_solver.expand_deflation_basis(U(Q));
        % Display of outputs.
        disp(['   lambda = ', num2str(lambda), ...
            ', d_lambda = ', num2str(d_lambda), ...
            ', omega = ', num2str(omega), ...
            ', d_omega = ', num2str(d_omega), ...
            ', d_lambda_diff_scaled = ', num2str(d_lambda/d_omega_test*(omega_hist(step) - omega_hist(1))), ...
            ', U_max = ', num2str(Umax_hist(step))]);
        
        if (d_lambda / d_omega_test) * (omega_hist(step) - omega_hist(1)) < d_lambda_diff_scaled_min
            disp('Minimal increase of lambda was achieved.');
            break;
        end
        
        % Update of d_lambda and d_omega.
        if d_omega_test > 1.5 * d_omega
            d_lambda = d_lambda / 2;
        end
        % if d_omega_test > 2 * d_omega
        %     d_lambda = d_lambda / 2;
        % end
        d_omega = d_omega_test;
    end % if flag
    
    % Stopping criteria.
    if d_lambda < d_lambda_min
        disp('Minimal increment of lambda was achieved.');
        break;
    end
    if step >= step_max
        disp('Maximal number of steps was achieved.');
        break;
    end
    
end % while lambda

% Clipping of the output arrays.
lambda_hist = lambda_hist(1:step);
omega_hist = omega_hist(1:step);
Umax_hist = Umax_hist(1:step);
work_hist = work_hist(1:step);

end % function
