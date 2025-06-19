function [U, lambda_hist, omega_hist, Umax_hist] = SSR_indirect_continuation(...
    lambda_init, d_lambda_init, d_lambda_min, d_lambda_diff_scaled_min, step_max, ...
    omega_max_stop, it_newt_max, it_damp_max, tol, r_min, ...
    K_elast, Q, f, constitutive_matrix_builder, linear_system_solver)
%--------------------------------------------------------------------------
% SSR_indirect_continuation performs an indirect continuation method to 
% estimate the safety factor lambda* and to construct the curve between the 
% strength reduction parameter lambda and the parameter omega representing 
% the work of external forces.
%
% This method incrementally increases omega and, at each step, solves the 
% nonlinear system:
%
%       F_lambda(U) = f   and   f' * U = omega,
%
% where F_lambda(U) is the internal force vector for a given lambda.
%
% The algorithm updates lambda and omega adaptively, storing the history of 
% computed values for further analysis.
%
% INPUTS:
%   lambda_init              - Initial value of the strength reduction parameter.
%   d_lambda_init            - Initial increment of lambda.
%   d_lambda_min             - Minimal increment of lambda.
%   d_lambda_diff_scaled_min - Minimal scaled increase in lambda (stopping criterion).
%   step_max                 - Maximum number of continuation steps.
%   omega_max_stop           - Maximum value of omega.
%   it_newt_max              - Maximum number of Newton iterations per continuation step.
%   it_damp_max              - Maximum number of damping iterations in Newton's method.
%   tol                      - Relative tolerance for Newton's solver.
%   r_min                    - Basic regularization parameter for the stiffness matrix.
%   K_elast                  - Elastic stiffness matrix.
%   Q                        - Logical array restricting nodes with Dirichlet boundary cond.
%   f                        - External load vector.
%   constitutive_matrix_builder - Object to build the constitutive force vector and tangent.
%   linear_system_solver     - Linear solver object (with deflation capability).
%
% OUTPUTS:
%   U           - Displacement field corresponding to lambda*.
%   lambda_hist - History of lambda values (for constructing the lambda-omega curve).
%   omega_hist  - History of omega values.
%   Umax_hist   - History of the maximum nodal displacement per step.
%
%--------------------------------------------------------------------------
%
% Initialization
%
lambda_hist = zeros(1, 1000);    % History of lambda values.
omega_hist = zeros(1, 1000);     % History of omega values.
Umax_hist  = zeros(1, 1000);      % History of maximum displacement values.

% First two steps of the continuation method.
[U_old, U, omega_old, omega, lambda_old, lambda] = CONTINUATION.init_phase_SSR_indirect_continuation(...
    lambda_init, d_lambda_init, d_lambda_min, ...
    it_newt_max, it_damp_max, tol, r_min, K_elast, Q, f, ...
    constitutive_matrix_builder, linear_system_solver.copy());

linear_system_solver.expand_deflation_basis(U_old(Q));

% Storage of the computed values.
omega_hist(1) = omega_old;
lambda_hist(1) = lambda_old;
omega_hist(2) = omega;
lambda_hist(2) = lambda;
Umax_hist(1) = max(sqrt(sum(U_old.^2, 1)));
Umax_hist(2) = max(sqrt(sum(U.^2, 1)));

% Other initial data.
d_omega = omega - omega_old;     % Current increment of omega.

%
% While cycle over the parameter omega.
%
step = 2;            % Number of continuation steps.
n_omega = 0;         % Number of reductions of omega.
n_omega_max = 5;     % Maximal number of reductions of omega.
while true
    fprintf('\n');
    fprintf(' Step = %d\n', step + 1);
    fprintf('\n');

    % Update of the parameter omega.
    omega_it = min(omega + d_omega, omega_max_stop);
    d_omega = omega_it - omega;

    % Initial estimate of the displacement field by linear extrapolation.
    U_ini = d_omega * (U - U_old) / (omega - omega_old) + U;

    % Solvers for the system F_lambda(U) = f, f'*U = omega_it.
    % a) Newton's method
    [U_it, lambda_it, flag, ~] = NEWTON.newton_ind_SSR(U_ini, omega_it, lambda, ...
        it_newt_max, it_damp_max, tol, r_min, K_elast, Q, f, ...
        constitutive_matrix_builder, linear_system_solver.copy());

    %       % b) bisection method w.r.t lambda
    %       [U_it, lambda_it, flag] = bisection_method(U_ini, omega_it, lambda, d_lambda, ...
    %           it_newt_max, it_damp_max, tol, r_min, r_damp, ...
    %           WEIGHT, B, K_elast, Q, f, c0, phi, psi, Davis_type, shear, bulk, lame);

    %       % c) tangent method w.r.t lambda
    %       [U_it, lambda_it, flag] = tangent_method(U_ini, omega_it, lambda, d_lambda, ...
    %           it_newt_max, it_damp_max, tol, r_min, r_damp, ...
    %           WEIGHT, B, K_elast, Q, f, c0, phi, psi, Davis_type, shear, bulk, lame);

    % Evaluation of the solver and update.
    if flag == 1  % Solver was not successful.
        d_omega = d_omega / 2;
        n_omega = n_omega + 1;
    else        % Solver was successful.
        % Update.
        step = step + 1;
        U_old = U;
        U = U_it;
        linear_system_solver.expand_deflation_basis(U(Q));
        omega_old = omega;
        omega = omega_it;
        d_lambda = lambda_it - lambda;
        lambda = lambda_it;
        n_omega = 0;
        % Storage.
        omega_hist(step) = omega;
        lambda_hist(step) = lambda;
        Umax_hist(step) = max(sqrt(sum(U.^2, 1)));
        % Display.
        disp(['   lambda = ', num2str(lambda), ...
            ', d_lambda = ', num2str(d_lambda), ...
            ', d_lambda_diff_scaled = ', num2str(d_lambda / d_omega * (omega_hist(step) - omega_hist(1))), ...
            ', omega = ', num2str(omega), ...
            ', d_omega = ', num2str(d_omega), ...
            ', U_max = ', num2str(Umax_hist(step))]);

        if (d_lambda / d_omega) * (omega_hist(step) - omega_hist(1)) < d_lambda_diff_scaled_min
            disp('Minimal increase of lambda was achieved.');
            break;
        end

        if omega >= omega_max_stop
            disp('Maximal omega was achieved.');
            break;
        end

        % Criterion for the change of d_omega.
        if (lambda_hist(step) - lambda_hist(step-1)) < 0.9 * (lambda_hist(step-1) - lambda_hist(step-2))
            d_omega = 2 * d_omega;
        end
    end % if flag

    % Stopping criteria for the indirect continuation method.
    if n_omega >= n_omega_max
        disp('It is impossible to increase omega.');
        break;
    end
    if step >= step_max
        disp('Maximal number of steps was achieved.');
        break;
    end

end % while

% Clipping of the output arrays.
lambda_hist = lambda_hist(1:step);
omega_hist  = omega_hist(1:step);
Umax_hist   = Umax_hist(1:step);

end % function
