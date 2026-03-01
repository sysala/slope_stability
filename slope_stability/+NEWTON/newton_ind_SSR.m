function [U_it, lambda_it, flag_N, it, history] = newton_ind_SSR(U_ini, omega, lambda_ini, ...
    it_newt_max, it_damp_max, tol, r_min, K_elast, Q, f, constitutive_matrix_builder, linear_solver)
%--------------------------------------------------------------------------
% newton_ind_SSR solves the system
%
%       F_lambda(U) = f   and   f' * U = omega,
%
% using a semismooth Newton method with regularization.
%
% OUTPUTS:
%   U_it      - Updated displacement field after convergence or termination.
%   lambda_it - Updated value of lambda.
%   flag_N    - Newton success/failure flag (0: success, 1: failure).
%   it        - Number of Newton iterations performed.
%   history   - Struct with fields:
%                 .residual  - History of relative residuals.
%                 .r         - History of regularization parameters.
%                 .alpha     - History of damping factors.
%                 .timing    - Accumulated wall-clock profiling.
%--------------------------------------------------------------------------

% Initialization.
n_n = size(U_ini, 2);
dim = size(U_ini, 1);
V = zeros(dim, n_n);
W = zeros(dim, n_n);
U_it = U_ini;
lambda_it = lambda_ini;
eps = tol / 1000;
norm_f = norm(f(Q));
flag_N = 0;
compute_diffs = 1;
r = r_min;
it = 0;

% Preallocate history arrays
residual_history = nan(1, it_newt_max);
r_history = nan(1, it_newt_max);
alpha_history = nan(1, it_newt_max);

% Timing accumulators for profiling
t_build_F_and_DS = 0;
t_build_F_eps = 0;
t_K_r_assembly = 0;
t_sparsersb_build = 0;
t_setup_preconditioner = 0;
t_A_orthogonalize = 0;
t_solve_W = 0;
t_solve_V = 0;
t_damping = 0;
t_expand_deflation_W = 0;
t_expand_deflation_V = 0;

% Sub-timing accumulators: build_F_and_DS_all breakdown
t_build_F_DS__reduction      = 0;
t_build_F_DS__stress_tangent = 0;
t_build_F_DS__build_F        = 0;
n_build_F_DS__calls          = 0;

% Sub-timing accumulators: damping_ALG5 breakdown
t_damping__build_F = 0;
t_damping__norm    = 0;
n_damping__iters   = 0;
n_damping__calls   = 0;

% Sub-timing accumulators: dfgmres_solver breakdown (V + W combined)
t_solve__precond     = 0;
t_solve__project     = 0;
t_solve__matvec      = 0;
t_solve__ortho       = 0;
t_solve__leastsq     = 0;
t_solve__init        = 0;
t_solve__reconstruct = 0;
n_solve__gmres_iters = 0;
n_solve__matvecs     = 0;
n_solve__prec_applies = 0;
n_solve__calls       = 0;

% One-time: restrict B to free DOFs and build maximal sparsity pattern.
if ~constitutive_matrix_builder.pattern_initialized
    constitutive_matrix_builder.init_K_r_pattern(Q);
end
n_Q = constitutive_matrix_builder.n_Q;

% Precompute upper-triangle mask for sparsersb "unique","sym" constructor.
% ref_I/ref_J are the same every call, so we mask once.
upper_mask = (constitutive_matrix_builder.ref_I <= constitutive_matrix_builder.ref_J);
K_I_sym = constitutive_matrix_builder.ref_I(upper_mask);
K_J_sym = constitutive_matrix_builder.ref_J(upper_mask);

% Flag for HYPRE IJV pattern reuse: first call does full setup,
% subsequent calls only update values (same sparsity pattern).
hypre_pattern_initialized = false;

% Semismooth Newton loop.
while true
    it = it + 1;
    t_step = tic;

    %% Compute constitutive response (F and DS, skip sparse K_tangent).
    if compute_diffs
        t_tmp = tic;
        F = constitutive_matrix_builder.build_F_and_DS_all(lambda_it, U_it);
        t_build_F_and_DS = t_build_F_and_DS + toc(t_tmp);
        n_build_F_DS__calls = n_build_F_DS__calls + 1;
        sub = constitutive_matrix_builder.last_build_F_DS_timing;
        if ~isempty(fieldnames(sub))
            t_build_F_DS__reduction      = t_build_F_DS__reduction      + sub.t_reduction;
            t_build_F_DS__stress_tangent = t_build_F_DS__stress_tangent + sub.t_stress_tangent;
            t_build_F_DS__build_F        = t_build_F_DS__build_F        + sub.t_build_F;
        end
        criterion = norm(F(Q) - f(Q));
        rel_resid = criterion / norm_f;

        residual_history(it) = rel_resid;

        if (rel_resid < tol) && (it > 1)
            fprintf('Newton method converges: iteration = %d, rel. resid = %e\n', it, rel_resid);
            break;
        end
    end

    r_history(it) = r;

    %% Build K_r(Q,Q) as [I, J, V] triplets, then sparsersb matrix.
    t_tmp = tic;
    [K_I, K_J, K_V] = constitutive_matrix_builder.build_K_r_QQ_vals(r);
    t_K_r_assembly = t_K_r_assembly + toc(t_tmp);

    t_tmp = tic;
    K_rQQ = sparsersb(K_I_sym, K_J_sym, K_V(upper_mask), n_Q, n_Q, "unique", "sym");
    t_sparsersb_build = t_sparsersb_build + toc(t_tmp);

    %% Estimate Newton increment.
    lambda_eps = lambda_it + eps;
    t_tmp = tic;
    F_eps = constitutive_matrix_builder.build_F_all(lambda_eps, U_it);
    t_build_F_eps = t_build_F_eps + toc(t_tmp);
    G = (F_eps - F) / eps;

    %% Setup/update preconditioner (HYPRE reuses sparsity pattern after first call).
    t_tmp = tic;
    if ~hypre_pattern_initialized
        linear_solver.setup_preconditioner_ijv(K_I, K_J, K_V, n_Q);
        hypre_pattern_initialized = true;
    else
        linear_solver.update_preconditioner_values(K_V);
    end
    t_setup_preconditioner = t_setup_preconditioner + toc(t_tmp);

    t_tmp = tic;
    linear_solver.A_orthogonalize(K_rQQ);
    t_A_orthogonalize = t_A_orthogonalize + toc(t_tmp);

    t_tmp = tic;
    W(Q) = linear_solver.solve(K_rQQ, -G(Q));
    t_solve_W = t_solve_W + toc(t_tmp);
    n_solve__calls = n_solve__calls + 1;
    st = linear_solver.last_solve_timing;
    if ~isempty(fieldnames(st))
        t_solve__precond     = t_solve__precond     + st.t_precond;
        t_solve__project     = t_solve__project     + st.t_project;
        t_solve__matvec      = t_solve__matvec      + st.t_matvec;
        t_solve__ortho       = t_solve__ortho       + st.t_ortho;
        t_solve__leastsq     = t_solve__leastsq     + st.t_leastsq;
        t_solve__init        = t_solve__init        + st.t_init;
        t_solve__reconstruct = t_solve__reconstruct + st.t_reconstruct;
        n_solve__gmres_iters = n_solve__gmres_iters + st.n_iters;
        n_solve__matvecs     = n_solve__matvecs     + st.n_matvecs;
        n_solve__prec_applies = n_solve__prec_applies + st.n_prec_applies;
    end
    t_tmp = tic;
    V(Q) = linear_solver.solve(K_rQQ, f(Q) - F(Q));
    t_solve_V = t_solve_V + toc(t_tmp);
    n_solve__calls = n_solve__calls + 1;
    st = linear_solver.last_solve_timing;
    if ~isempty(fieldnames(st))
        t_solve__precond     = t_solve__precond     + st.t_precond;
        t_solve__project     = t_solve__project     + st.t_project;
        t_solve__matvec      = t_solve__matvec      + st.t_matvec;
        t_solve__ortho       = t_solve__ortho       + st.t_ortho;
        t_solve__leastsq     = t_solve__leastsq     + st.t_leastsq;
        t_solve__init        = t_solve__init        + st.t_init;
        t_solve__reconstruct = t_solve__reconstruct + st.t_reconstruct;
        n_solve__gmres_iters = n_solve__gmres_iters + st.n_iters;
        n_solve__matvecs     = n_solve__matvecs     + st.n_matvecs;
        n_solve__prec_applies = n_solve__prec_applies + st.n_prec_applies;
    end

    d_l = - (f(Q)' * V(Q)) / (f(Q)' * W(Q));
    d_U = V + d_l * W;

    %% Line search damping.
    t_tmp = tic;
    [alpha, damp_timing] = NEWTON.damping_ALG5(it_damp_max, U_it, lambda_it, d_U, d_l, f, criterion, Q, constitutive_matrix_builder);
    t_damping = t_damping + toc(t_tmp);
    t_damping__build_F = t_damping__build_F + damp_timing.t_build_F;
    t_damping__norm    = t_damping__norm    + damp_timing.t_norm;
    n_damping__iters   = n_damping__iters   + damp_timing.n_iters;
    n_damping__calls   = n_damping__calls   + 1;
    alpha_history(it) = alpha;

    %% Regularization adjustments.
    compute_diffs = 1;
    if alpha < 1e-1
        if alpha == 0
            compute_diffs = 0;
            r = r * 2;
        else
            r = r * 2^(1/4);
        end
    else
        t_tmp = tic;
        linear_solver.expand_deflation_basis(W(Q));
        t_expand_deflation_W = t_expand_deflation_W + toc(t_tmp);
        t_tmp = tic;
        linear_solver.expand_deflation_basis(V(Q));
        t_expand_deflation_V = t_expand_deflation_V + toc(t_tmp);
        if alpha > 0.5
            r = max(r / sqrt(2), r_min);
        end
    end

    if (alpha == 0) && (r > 1)
        fprintf('\nNewton solver does not converge: iteration = %d, stopping criterion=%e\n', it, rel_resid);
        if rel_resid > 10 * tol
            flag_N = 1;
        end
        break;
    end

    %% Update variables.
    U_it = U_it + alpha * d_U;
    U_it = omega * U_it / (f(Q)' * U_it(Q));  % Enforce constraint
    lambda_it = lambda_it + alpha * d_l;

    fprintf('  newton_ind_SSR it=%d  resid=%.2e  alpha=%.2f  r=%.3g  step_time=%.2f s\n', ...
        it, rel_resid, alpha, r, toc(t_step));

    %% Check maximum iterations.
    if it == it_newt_max
        fprintf('Newton solver does not converge: iteration = %d, stopping criterion=%e\n', it, rel_resid);
        if rel_resid > 10 * tol
            flag_N = 1;
        end
        break;
    end
end % while

% Trim unused history entries
history.residual = residual_history(1:it);
history.r = r_history(1:it);
history.alpha = alpha_history(1:it);

% Store timing profile
history.timing.build_F_and_DS      = t_build_F_and_DS;
history.timing.solve_V             = t_solve_V;
history.timing.A_orthogonalize     = t_A_orthogonalize;
history.timing.damping             = t_damping;
history.timing.setup_preconditioner = t_setup_preconditioner;
history.timing.build_F_eps         = t_build_F_eps;
history.timing.solve_W             = t_solve_W;
history.timing.K_r_assembly        = t_K_r_assembly;
history.timing.sparsersb_build     = t_sparsersb_build;
history.timing.expand_deflation_W  = t_expand_deflation_W;
history.timing.expand_deflation_V  = t_expand_deflation_V;

% Sub-profiling: build_F_and_DS_all breakdown
history.timing.build_F_DS__reduction      = t_build_F_DS__reduction;
history.timing.build_F_DS__stress_tangent = t_build_F_DS__stress_tangent;
history.timing.build_F_DS__build_F        = t_build_F_DS__build_F;
history.timing.build_F_DS__n_calls        = n_build_F_DS__calls;

% Sub-profiling: damping_ALG5 breakdown
history.timing.damping__build_F = t_damping__build_F;
history.timing.damping__norm    = t_damping__norm;
history.timing.damping__n_iters = n_damping__iters;
history.timing.damping__n_calls = n_damping__calls;

% Sub-profiling: dfgmres_solver breakdown (V + W combined)
history.timing.solve__precond      = t_solve__precond;
history.timing.solve__project      = t_solve__project;
history.timing.solve__matvec       = t_solve__matvec;
history.timing.solve__ortho        = t_solve__ortho;
history.timing.solve__leastsq      = t_solve__leastsq;
history.timing.solve__init         = t_solve__init;
history.timing.solve__reconstruct  = t_solve__reconstruct;
history.timing.solve__n_gmres_iters = n_solve__gmres_iters;
history.timing.solve__n_matvecs    = n_solve__matvecs;
history.timing.solve__n_prec_applies = n_solve__prec_applies;
history.timing.solve__n_calls      = n_solve__calls;

end % function
