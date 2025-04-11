function [U2, omega1, omega2, lambda2, U1] = init_phase_SSR_direct_continuation(...
                                 lambda1, d_lambda1, d_lambda_min, ...
                                 it_newt_max, it_damp_max, tol, eps, r_min, K_elast, Q, f, ...
                                 constitutive_matrix_builder, linear_system_solver)
%--------------------------------------------------------------------------
% init_phase_SSR_direct_continuation performs the initial phase for the 
% continuation method (ALG2). The aim is to compute the remaining input data 
% required by the continuation algorithm.
%
% In particular, given an initial strength reduction parameter lambda1 and an 
% initial increment d_lambda1, this function computes:
%
%   U2      - The displacement field corresponding to a new parameter value lambda2.
%   omega1  - The value of the control parameter omega corresponding to lambda1.
%   omega2  - The value of omega corresponding to lambda2.
%   lambda2 - The second value of the parameter lambda for which the Newton method converges.
%
% The procedure is as follows:
%   1. For the first step, the function computes U1 and omega1 by solving 
%      F_lambda(U) = f (with lambda = lambda1) starting from an initial guess U_ini.
%   2. In the second step, the parameter lambda is updated to lambda1 + d_lambda.
%      The function then attempts to compute a new solution U and omega.
%      If the solver fails (or if the computed omega does not exceed omega1),
%      d_lambda is halved until a successful solution is obtained.
%
% INPUTS:
%   lambda1                 - Initial value of the parameter lambda.
%   d_lambda1               - Initial increment of lambda.
%   d_lambda_min            - Minimal allowed increment of lambda.
%   it_newt_max             - Maximum number of Newton iterations.
%   it_damp_max             - Maximum number of damping iterations in Newton's method.
%   tol                     - Relative tolerance for Newton's solvers.
%   eps                     - Precision for numerical derivatives.
%   r_min                   - Basic regularization parameter for the stiffness matrix.
%   K_elast                 - Elastic stiffness matrix.
%   Q                       - Logical array restricting nodes with Dirichlet boundary conditions.
%   f                       - External load vector.
%   constitutive_matrix_builder - Object to build the constitutive force vector and tangent operator.
%   linear_system_solver    - Linear solver object (with deflation capabilities).
%
% OUTPUTS:
%   U2      - Displacement field corresponding to lambda2.
%   omega1  - Value of the control parameter omega corresponding to lambda1.
%   omega2  - Value of the control parameter omega corresponding to lambda2.
%   lambda2 - The updated value of the parameter lambda.
%
%--------------------------------------------------------------------------
%
% Initialization.
n_n = size(Q, 2);                % Number of nodes.
dim = size(Q, 1);                % Spatial dimension.
U_ini = zeros(dim, n_n);         % Initial displacement field.

%
% First step of the continuation method.
%
fprintf('\n'); 
fprintf(' Step = %d  ', 1); 
fprintf('\n'); 

% Compute U and omega for the initial lambda (lambda1).
[U1, omega1, flag] = CONTINUATION.omega_SSR_direct_continuation(...
    lambda1, U_ini, eps, 1000, ...
    it_newt_max, it_damp_max, tol, r_min, K_elast, Q, f, ...
    constitutive_matrix_builder, linear_system_solver.copy());

% Test convergence: if the solver fails, the initial lambda is too large.
if flag == 1
    error('Initial choice of lambda seems to be too large.');
end

disp(['   lambda = ', num2str(lambda1), ...
      ', d_lambda = ', num2str(d_lambda1), ...
      ', omega = ', num2str(omega1)]);

%
% Second step of the continuation method.
%
fprintf('\n'); 
fprintf(' Step = %d  ', 2); 
fprintf('\n'); 

d_lambda = d_lambda1;  % Set the initial increment of lambda.
while true                   
                   
    % Update parameter lambda.
    lambda_it = lambda1 + d_lambda;
                
    % Expand deflation basis with U1 restricted to Q.
    linear_system_solver.expand_deflation_basis(U1(Q));
    
    % Compute U and omega for lambda_it.
    [U_it, omega_it, flag] = CONTINUATION.omega_SSR_direct_continuation(...
        lambda_it, U1, eps, d_lambda, ...
        it_newt_max, it_damp_max, tol, r_min, K_elast, Q, f, ...
        constitutive_matrix_builder, linear_system_solver.copy());
        
    % Test convergence.
    if (flag == 1) || (omega_it <= omega1)
        % If the solver was not successful or omega did not increase, decrease d_lambda.
        d_lambda = d_lambda / 2;
    else
        % If the solver was successful, accept the new parameters.
        U2 = U_it;
        omega2 = omega_it;
        lambda2 = lambda_it;  
        break;
    end

    % Stopping criteria: if d_lambda is too small.
    if d_lambda < d_lambda_min
        error('It seems that FoS is equal to lambda_init.');
    end
end

disp(['   lambda = ', num2str(lambda2), ...
      ', d_lambda = ', num2str(lambda2 - lambda1), ...
      ', omega = ', num2str(omega1),...
      ', d_omega = ', num2str(omega2 - omega1)]);

end % function
