function [U_it, flag_N, it, history] = newton(U_ini, tol, it_newt_max, it_damp_max, r_min, ...
    K_elast, Q, f, constitutive_matrix_builder, linear_system_solver)
%--------------------------------------------------------------------------
% newton applies a semismooth Newton method with damping and regularization
% to solve the nonlinear system:
%
%       F(U) = f,
%
% where F(U) is the internal force vector computed from the displacement field U.
%
% The method updates U iteratively until the relative residual
% norm(F(Q)-f(Q))/norm(f(Q)) is below the prescribed tolerance tol.
%
% INPUT ARGUMENTS:
%   U_ini               - Initial displacement field, size (dim, n_n).
%   tol                 - Relative tolerance for Newton's convergence.
%   it_newt_max         - Maximum number of Newton iterations.
%   it_damp_max         - Maximum number of damping iterations.
%   r_min               - Base regularization parameter.
%   K_elast             - (unused, kept for API compatibility; elastic K
%                          is precomputed inside constitutive_matrix_builder).
%   Q                   - Logical free-DOF mask (dim x n_n).
%   f                   - External load vector.
%   constitutive_matrix_builder - CONSTITUTIVE object (pattern must be
%                         initialised via init_K_r_pattern before calling).
%   linear_system_solver- DFGMRES solver object with IJV support.
%
% OUTPUTS:
%   U_it    - Computed displacement field approximating the solution.
%   flag_N  - Newton solver flag: 0 if converged, 1 if not converged.
%   it      - Number of Newton iterations performed.
%   history - Struct with timing profile (same fields as newton_ind_SSR).
%--------------------------------------------------------------------------
%
% Initialization.
n_n = size(U_ini, 2);         % Number of nodes.
dim = size(U_ini, 1);         % Spatial dimension.
dU = zeros(dim, n_n);         % Newton increment.
U_it = U_ini;                 % Current displacement.
it = 0;                       % Iteration counter.
flag_N = 0;                   % Convergence flag.
r = r_min;                    % Regularization parameter.
compute_diffs = 1;            % Flag: must recompute stress/tangent?

% Ensure pattern is initialised.
if ~constitutive_matrix_builder.pattern_initialized
    constitutive_matrix_builder.init_K_r_pattern(Q);
end
n_Q = constitutive_matrix_builder.n_Q;

% Timing accumulators.
t_build_F_and_DS = 0;
t_K_r_assembly   = 0;
t_sparsersb_build = 0;
t_setup_preconditioner = 0;
t_A_orthogonalize = 0;
t_solve = 0;
t_damping = 0;
t_expand_deflation = 0;

% Newton's iteration loop.
while true
    it = it + 1;
    t_step = tic;

    % Compute internal force and tangent moduli.
    if compute_diffs
        t_tmp = tic;
        F = constitutive_matrix_builder.build_F_and_DS_reduced(U_it);
        t_build_F_and_DS = t_build_F_and_DS + toc(t_tmp);
    end

    % Compute relative residual criterion.
    criterion = norm(F(Q) - f(Q)) / norm(f(Q));
    if (criterion < tol)
        fprintf('Newton method converges: iteration = %d, stopping criterion = %e\n', it, criterion);
        break;
    end

    % Build K_r(Q,Q) as [I, J, V] triplets  →  sparsersb.
    t_tmp = tic;
    [K_I, K_J, K_V] = constitutive_matrix_builder.build_K_r_QQ_vals(r);
    t_K_r_assembly = t_K_r_assembly + toc(t_tmp);

    t_tmp = tic;
    K_rQQ = sparsersb(K_I, K_J, K_V, n_Q, n_Q);
    t_sparsersb_build = t_sparsersb_build + toc(t_tmp);

    % Setup preconditioner  (HYPRE reuses pattern automatically on same instance).
    t_tmp = tic;
    linear_system_solver.setup_preconditioner_ijv(K_I, K_J, K_V, n_Q);
    t_setup_preconditioner = t_setup_preconditioner + toc(t_tmp);

    t_tmp = tic;
    linear_system_solver.A_orthogonalize(K_rQQ);
    t_A_orthogonalize = t_A_orthogonalize + toc(t_tmp);

    % Solve for the Newton increment.
    t_tmp = tic;
    dU(Q) = linear_system_solver.solve(K_rQQ, f(Q) - F(Q));
    t_solve = t_solve + toc(t_tmp);

    % Determine damping factor via line search.
    t_tmp = tic;
    alpha = NEWTON.damping(it_damp_max, U_it, dU, F, f, constitutive_matrix_builder);
    t_damping = t_damping + toc(t_tmp);

    % Regularization adjustments based on alpha.
    compute_diffs = 1;
    if alpha < 1e-1
        if alpha == 0
            compute_diffs = 0;  % Skip recomputation if step is rejected.
            r = r * 2;
        else
            r = r * 2^(1/4);
        end
    else
        % Expand deflation basis with the computed increment.
        t_tmp = tic;
        linear_system_solver.expand_deflation_basis(dU(Q));
        t_expand_deflation = t_expand_deflation + toc(t_tmp);
        if alpha > 0.5
            r = max(r / sqrt(2), r_min);
        end
    end

    if (alpha == 0) && (r > 1)
        fprintf('\nNewton solver does not converge: stopping criterion = %e\n', criterion);
        flag_N = 1;
        break;
    end

    % Update displacement.
    U_it = U_it + alpha * dU;

    fprintf('  newton it=%d  resid=%.2e  alpha=%.2f  r=%.3g  step_time=%.2f s\n', ...
        it, criterion, alpha, r, toc(t_step));

    % Check for divergence or maximum iterations.
    if isnan(criterion) || (it == it_newt_max)
        flag_N = 1;
        fprintf('Newton solver does not converge: stopping criterion = %e\n', criterion);
        break;
    end

end % while

% Build timing profile.
history.timing.build_F_and_DS       = t_build_F_and_DS;
history.timing.K_r_assembly         = t_K_r_assembly;
history.timing.sparsersb_build      = t_sparsersb_build;
history.timing.setup_preconditioner = t_setup_preconditioner;
history.timing.A_orthogonalize      = t_A_orthogonalize;
history.timing.solve                = t_solve;
history.timing.damping              = t_damping;
history.timing.expand_deflation     = t_expand_deflation;

end % function
