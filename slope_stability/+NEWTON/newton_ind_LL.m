function [U_it, t_it, flag_N, it] = newton_ind_LL( ...
    U_ini, t_ini, omega, it_newt_max, it_damp_max, tol, r_min, K_elast, Q, f, ...
    constitutive_matrix_builder, linear_solver)
%--------------------------------------------------------------------------
% newton_ind_LL applies a semismooth Newton method to solve:
%
%       F(U) = t * f,    and    f' * U = omega,
%
% where F(U) is the internal force vector depending on the displacement U.
%
% The method adaptively updates both the displacement field U and the load
% factor t, employing a line search damping strategy and regularization
% based on a mixture of the elastic stiffness and the consistent tangent
% operator.
%
% INPUT ARGUMENTS:
%   U_ini               - Initial displacement field.
%   t_ini               - Initial load factor.
%   omega               - Prescribed control parameter.
%   it_newt_max         - Maximum number of Newton iterations.
%   it_damp_max         - Maximum number of damping (line search) steps.
%   tol                 - Convergence tolerance (relative).
%   r_min               - Base regularization factor for the stiffness matrix.
%   K_elast             - Elastic stiffness matrix.
%   Q                   - Logical array (mask) for restricted Dirichlet dofs.
%   f                   - External load vector.
%   constitutive_matrix_builder - Object to build constitutive matrices.
%   linear_solver       - Linear solver object with deflation capabilities.
%
% OUTPUT:
%   U_it  - Updated displacement field upon convergence or termination.
%   t_it  - Updated load factor.
%   flag_N- Newton success/failure flag (0: success, 1: failure).
%   it    - Number of Newton iterations performed.
%
%--------------------------------------------------------------------------
%
% Basic size and geometry definitions.
n_n = size(U_ini, 2);    % Number of nodes.
dim = size(U_ini, 1);    % Dimension (2D or 3D).

% Initialize auxiliary vectors for iterative updates.
V = zeros(dim, n_n);
W = zeros(dim, n_n);

% Initialize iteration variables.
U_it  = U_ini;         % Current displacement.
t_it  = t_ini;         % Current load factor.
flag_N = 0;            % Success flag: 0 if converged, 1 if failure.
it    = 0;             % Newton iteration counter.
norm_f = norm(f(Q));   % Norm of the restricted external load.
r = r_min;             % Regularization parameter.

compute_diffs = 1;     % Flag to decide whether to recompute differences.

% Start semismooth Newton iteration.
while true
    it = it + 1;
    
    % Compute internal force vector and tangent operator.
    if compute_diffs
        [F_int, K_tangent] = constitutive_matrix_builder.build_F_K_tangent_reduced(U_it);
    else
        F_int = constitutive_matrix_builder.build_F_reduced(U_it);
    end

    % Compute convergence criterion.
    criterion = norm(F_int(Q) - t_it * f(Q));
    if (criterion < tol * norm_f) && (it > 1)
        fprintf('Newton converged in %d iterations. Rel. resid = %e\n', it, criterion / norm_f);
        break;
    end

    % Form the regularized stiffness matrix K_r as a convex combination of K_elast and K_tangent.
    K_r = r * K_elast + (1 - r) * K_tangent;
    
    % Set up and orthogonalize the preconditioner for the restricted system.
    linear_solver.setup_preconditioner(K_r(Q, Q));
    linear_solver.A_orthogonalize(K_r(Q, Q));
    
    % Solve for W(Q): used for load increment computation.
    W(Q) = linear_solver.solve(K_r(Q, Q), f(Q));
    % Solve for V(Q): represents the displacement correction.
    V(Q) = linear_solver.solve(K_r(Q, Q), t_it * f(Q) - F_int(Q));
    
    % Compute the load increment d_t.
    d_t = - (f(Q)' * V(Q)) / (f(Q)' * W(Q));
    % Compute the displacement increment dU.
    dU  = V + d_t * W;
    
    % Perform a line search damping to obtain a suitable step length alpha.
    alpha = NEWTON.damping(it_damp_max, U_it, dU, F_int, F_int*0, constitutive_matrix_builder);
    
    % Adjust regularization based on the damping result.
    compute_diffs = 1;
    if alpha < 1e-1
        if alpha == 0
            compute_diffs = 0;  % Skip recomputing K_tangent if step rejected.
            r = r * 2;
        else
            r = r * 2^(1/4);
        end
    else
        % Expand deflation bases with current correction vectors.
        linear_solver.expand_deflation_basis(W(Q));
        linear_solver.expand_deflation_basis(V(Q));
        
        if alpha > 0.5
            r = max(r / sqrt(2), r_min);
        end
    end
    if alpha == 0 && r > 1
        fprintf('\nNewton solver does not converge: stopping criterion=%e\n', criterion / norm_f);
        flag_N = 1;
        break;
    end

    % Update displacement and load factor.
    U_it = U_it + alpha * dU;
    % Enforce the constraint f' * U = omega.
    U_it = omega * U_it / (f(Q)' * U_it(Q));
    t_it = t_it + d_t;
    
    % Check if maximum iterations reached.
    if it == it_newt_max
        fprintf('Newton solver does not converge: stopping criterion=%e\n', criterion / norm_f);
        flag_N = 1;
        break;
    end
end % while

end % function
