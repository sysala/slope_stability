function [U, omega, flag] = omega_SSR_direct_continuation(lambda_ini, U_ini, eps, d_lambda, ...
    it_newt_max, it_damp_max, tol, r_min, K_elast, Q, f, ...
    constitutive_matrix_builder, linear_system_solver)
%--------------------------------------------------------------------------
% omega_SSR_direct_continuation computes the displacement field U and the 
% parameter omega for a given lambda, according to ALG2.
%
% This function solves the nonlinear system
%
%       F_lambda(U) = f,    and    omega = -G'(lambda),
%
% where F_lambda(U) denotes the internal force vector associated with the 
% strength reduction parameter lambda, and omega represents the derivative 
% (with a minus sign) of the cost function G with respect to lambda.
%
% The method first computes the minimizer U of the cost function J_lambda(U)
% and then computes an approximate numerical derivative in lambda by solving
% the problem at lambda - eps.
%
% INPUTS:
%   lambda_ini              - Current value of the parameter lambda.
%   U_ini                   - Initial guess for the displacement vector U.
%   eps                     - Precision for the numerical derivative (perturbation in lambda).
%   d_lambda                - Current increment in lambda.
%   it_newt_max             - Maximum number of Newton iterations.
%   it_damp_max             - Maximum number of damping iterations in Newton's method.
%   tol                     - Relative tolerance for Newton's solver.
%   r_min                   - Basic regularization parameter for the stiffness matrix.
%   K_elast                 - Elastic stiffness matrix.
%   Q                       - Logical array restricting nodes with Dirichlet boundary conditions.
%   f                       - External load vector.
%   constitutive_matrix_builder - Object to build the constitutive force vector and tangent;
%                                 it also performs material reduction.
%   linear_system_solver    - Linear solver object (with deflation capability).
%
% OUTPUTS:
%   U      - Computed displacement field (minimizer of J_lambda).
%   omega  - Computed value of the parameter omega = -G'(lambda).
%   flag   - Logical flag indicating Newton's convergence (0: success, 1: failure).
%
% NOTES:
%   - If the Newton solver does not converge (flag==1), omega remains at its
%     initial unrealistic value (0).
%
%--------------------------------------------------------------------------
%
% Initialization
%
omega = 0;  % Set an unrealistic value for omega in case the Newton solver fails.

% Compute the minimizer U of the cost function J_lambda and the corresponding flag.
constitutive_matrix_builder.reduction(lambda_ini);
[U, flag] = NEWTON.newton(U_ini, tol, it_newt_max, it_damp_max, r_min, K_elast, Q, f, ...
    constitutive_matrix_builder, linear_system_solver);

if flag == 1
    return
end

% Compute the integrated potential Psi and cost function J for the current U.
Psi_integrated = constitutive_matrix_builder.potential(U);
J = Psi_integrated - f(Q)' * U(Q);

% Compute the minimizer U_eps for the perturbed parameter (lambda - eps).
beta = min(1, eps / d_lambda);
U_beta = beta * U_ini + (1 - beta) * U;

constitutive_matrix_builder.reduction(lambda_ini - eps);
[U_eps, flag] = NEWTON.newton(U_beta, tol, it_newt_max, it_damp_max, r_min, K_elast, Q, f, ...
    constitutive_matrix_builder, linear_system_solver);

if flag == 1
    return
end

% Compute the cost function J_eps for U_eps.
Psi_integrated = constitutive_matrix_builder.potential(U_eps);
J_eps = Psi_integrated - f(Q)' * U_eps(Q);

% Compute omega as the numerical derivative of J with respect to lambda.
omega = (J_eps - J) / eps;

end % function
