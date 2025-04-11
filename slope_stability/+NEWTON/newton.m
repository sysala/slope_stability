function [U_it, flag_N, it] = newton(U_ini, tol, it_newt_max, it_damp_max, r_min, ...
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
%   U_ini               - Initial displacement field, size (dim, n_n)
%                         (dim: spatial dimension, n_n: number of nodes).
%   tol                 - Relative tolerance for Newton's convergence.
%   it_newt_max         - Maximum number of Newton iterations.
%   it_damp_max         - Maximum number of iterations for the line search
%                         damping procedure.
%   r_min               - Base regularization parameter for the tangent stiffness matrix.
%   K_elast             - Elastic stiffness matrix.
%   Q                   - Logical array (mask) for nodes with Dirichlet boundary conditions.
%   f                   - External load vector.
%   constitutive_matrix_builder - Object for building the constitutive force vector
%                                 and consistent tangent operator.
%   linear_system_solver- Linear solver object with preconditioning and deflation features.
%
% OUTPUTS:
%   U_it    - Computed displacement field approximating the solution.
%   flag_N  - Newton solver flag: 0 if converged, 1 if not converged.
%   it      - Number of Newton iterations performed.
%
%--------------------------------------------------------------------------
%
% Initialization.
n_n = size(U_ini, 2);         % Number of nodes.
dim = size(U_ini, 1);         % Spatial dimension (2D or 3D).
dU = zeros(dim, n_n);         % Newton increment.
U_it = U_ini;                 % Current displacement.
it = 0;                       % Iteration counter.
flag_N = 0;                   % Convergence flag.
r = r_min;                    % Regularization parameter.
compute_diffs = 1;            % Flag to indicate if K_tangent must be recomputed.

% Newton's iteration loop.
while true
    it = it + 1;
    
    % Compute internal force and tangent operator.
    if compute_diffs
        [F, K_tangent] = constitutive_matrix_builder.build_F_K_tangent_reduced(U_it);
    else
        F = constitutive_matrix_builder.build_F_reduced(U_it);
    end
    
    % Compute relative residual criterion.
    criterion = norm(F(Q) - f(Q)) / norm(f(Q));
    if (criterion < tol)
        fprintf('Newton method converges: iteration = %d, stopping criterion = %e\n', it, criterion);
        break;
    end
    
    % Form the regularized tangent stiffness matrix.
    K_r = r * K_elast + (1 - r) * K_tangent;
    
    % Solve for the Newton increment.
    linear_system_solver.setup_preconditioner(K_r(Q, Q));
    linear_system_solver.A_orthogonalize(K_r(Q, Q));
    dU(Q) = linear_system_solver.solve(K_r(Q, Q), f(Q) - F(Q));
    
    % Determine damping factor via line search.
    alpha = NEWTON.damping(it_damp_max, U_it, dU, F, f, constitutive_matrix_builder);
    
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
        linear_system_solver.expand_deflation_basis(dU(Q));
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
    
    % Check for divergence or maximum iterations.
    if isnan(criterion) || (it == it_newt_max)
        flag_N = 1;
        fprintf('Newton solver does not converge: stopping criterion = %e\n', criterion);
        break;
    end

end % while

end % function
