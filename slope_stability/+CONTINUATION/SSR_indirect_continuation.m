function [U, lambda_hist, omega_hist, Umax_hist, stats] = SSR_indirect_continuation(...
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
% where F_lambda(U) is the internal force vector for a given lambda. This
% system is primarily solved by the Newton method. Alternative methods 
% (e.g., bisection or tangent methods w.r.t. lambda) are provided as 
% commented-out code.
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
t_total = tic;

stats = struct();
stats.init_newton_iterations = [];
stats.init_linear_iterations = 0;
stats.init_linear_solve_time = 0;
stats.init_linear_preconditioner_time = 0;
stats.init_linear_orthogonalization_time = 0;
stats.attempt_step = [];
stats.attempt_success = [];
stats.attempt_wall_time = [];
stats.attempt_newton_iterations = [];
stats.attempt_newton_flag = [];
stats.attempt_newton_relres_end = [];
stats.attempt_linear_iterations = [];
stats.attempt_linear_solve_time = [];
stats.attempt_linear_preconditioner_time = [];
stats.attempt_linear_orthogonalization_time = [];
stats.attempt_omega_target = [];
stats.attempt_lambda_before = [];
stats.attempt_lambda_after = [];
stats.step_index = [];
stats.step_attempt_count = [];
stats.step_wall_time = [];
stats.step_newton_iterations = [];
stats.step_newton_iterations_total = [];
stats.step_newton_relres_end = [];
stats.step_linear_iterations = [];
stats.step_linear_solve_time = [];
stats.step_linear_preconditioner_time = [];
stats.step_linear_orthogonalization_time = [];
stats.step_lambda = [];
stats.step_omega = [];
stats.total_wall_time = 0;

lambda_hist = zeros(1, 1000);    % History of lambda values.
omega_hist = zeros(1, 1000);     % History of omega values.
Umax_hist  = zeros(1, 1000);      % History of maximum displacement values.

% First two steps of the continuation method.
snap_init_0 = local_collector_snapshot(linear_system_solver);
[U_old, U, omega_old, omega, lambda_old, lambda, init_newton_its] = CONTINUATION.init_phase_SSR_indirect_continuation(...
    lambda_init, d_lambda_init, d_lambda_min, ...
    it_newt_max, it_damp_max, tol, r_min, K_elast, Q, f, ...
    constitutive_matrix_builder, linear_system_solver.copy());
snap_init_1 = local_collector_snapshot(linear_system_solver);
delta_init = local_collector_delta(snap_init_0, snap_init_1);
stats.init_newton_iterations = init_newton_its;
stats.init_linear_iterations = delta_init.iterations;
stats.init_linear_solve_time = delta_init.solve_time;
stats.init_linear_preconditioner_time = delta_init.preconditioner_time;
stats.init_linear_orthogonalization_time = delta_init.orthogonalization_time;

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

% Check whether omega_max_stop is not too low
if omega_max_stop<omega + d_omega
    error('Too small value of omega_max_stop. It is necessary to increase this input quantinty.')
end

%
% While cycle over the parameter omega.
%
step = 2;            % Number of continuation steps.
n_omega = 0;         % Number of reductions of omega.
n_omega_max = 5;     % Maximal number of reductions of omega.
step_wall_accum = 0;
step_lin_it_accum = 0;
step_lin_solve_accum = 0;
step_lin_prec_accum = 0;
step_lin_orth_accum = 0;
step_newton_it_accum = 0;
step_attempt_count = 0;
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
    t_attempt = tic;
    snap_0 = local_collector_snapshot(linear_system_solver);
    [U_it, lambda_it, flag, it_newt, history] = NEWTON.newton_ind_SSR(U_ini, omega_it, lambda, ...
        it_newt_max, it_damp_max, tol, r_min, K_elast, Q, f, ...
        constitutive_matrix_builder, linear_system_solver.copy());
    snap_1 = local_collector_snapshot(linear_system_solver);
    delta_attempt = local_collector_delta(snap_0, snap_1);
    attempt_wall = toc(t_attempt);

    attempt_relres = NaN;
    if isfield(history, 'residual') && ~isempty(history.residual)
        attempt_relres = history.residual(end);
    end

    stats.attempt_step(end + 1, 1) = step + 1;
    stats.attempt_success(end + 1, 1) = (flag == 0);
    stats.attempt_wall_time(end + 1, 1) = attempt_wall;
    stats.attempt_newton_iterations(end + 1, 1) = it_newt;
    stats.attempt_newton_flag(end + 1, 1) = flag;
    stats.attempt_newton_relres_end(end + 1, 1) = attempt_relres;
    stats.attempt_linear_iterations(end + 1, 1) = delta_attempt.iterations;
    stats.attempt_linear_solve_time(end + 1, 1) = delta_attempt.solve_time;
    stats.attempt_linear_preconditioner_time(end + 1, 1) = delta_attempt.preconditioner_time;
    stats.attempt_linear_orthogonalization_time(end + 1, 1) = delta_attempt.orthogonalization_time;
    stats.attempt_omega_target(end + 1, 1) = omega_it;
    stats.attempt_lambda_before(end + 1, 1) = lambda;
    if flag == 0
        stats.attempt_lambda_after(end + 1, 1) = lambda_it;
    else
        stats.attempt_lambda_after(end + 1, 1) = NaN;
    end

    step_wall_accum = step_wall_accum + attempt_wall;
    step_lin_it_accum = step_lin_it_accum + delta_attempt.iterations;
    step_lin_solve_accum = step_lin_solve_accum + delta_attempt.solve_time;
    step_lin_prec_accum = step_lin_prec_accum + delta_attempt.preconditioner_time;
    step_lin_orth_accum = step_lin_orth_accum + delta_attempt.orthogonalization_time;
    step_newton_it_accum = step_newton_it_accum + it_newt;
    step_attempt_count = step_attempt_count + 1;

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

        stats.step_index(end + 1, 1) = step;
        stats.step_attempt_count(end + 1, 1) = step_attempt_count;
        stats.step_wall_time(end + 1, 1) = step_wall_accum;
        stats.step_newton_iterations(end + 1, 1) = it_newt;
        stats.step_newton_iterations_total(end + 1, 1) = step_newton_it_accum;
        stats.step_newton_relres_end(end + 1, 1) = attempt_relres;
        stats.step_linear_iterations(end + 1, 1) = step_lin_it_accum;
        stats.step_linear_solve_time(end + 1, 1) = step_lin_solve_accum;
        stats.step_linear_preconditioner_time(end + 1, 1) = step_lin_prec_accum;
        stats.step_linear_orthogonalization_time(end + 1, 1) = step_lin_orth_accum;
        stats.step_lambda(end + 1, 1) = lambda;
        stats.step_omega(end + 1, 1) = omega;

        step_wall_accum = 0;
        step_lin_it_accum = 0;
        step_lin_solve_accum = 0;
        step_lin_prec_accum = 0;
        step_lin_orth_accum = 0;
        step_newton_it_accum = 0;
        step_attempt_count = 0;

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
stats.total_wall_time = toc(t_total);

end % function

function snap = local_collector_snapshot(linear_system_solver)
collector = linear_system_solver.iteration_collector;
snap.iterations = collector.get_total_iterations();
snap.solve_time = collector.get_total_solve_time();
snap.preconditioner_time = collector.get_total_preconditioner_time();
snap.orthogonalization_time = collector.get_total_orthogonalization_time();
end

function delta = local_collector_delta(snap0, snap1)
delta.iterations = snap1.iterations - snap0.iterations;
delta.solve_time = snap1.solve_time - snap0.solve_time;
delta.preconditioner_time = snap1.preconditioner_time - snap0.preconditioner_time;
delta.orthogonalization_time = snap1.orthogonalization_time - snap0.orthogonalization_time;
end
