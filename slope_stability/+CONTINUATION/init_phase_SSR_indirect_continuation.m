function [U1, U2, omega1, omega2, lambda1, lambda2, all_newton_its] = init_phase_SSR_indirect_continuation(...
                                 lambda_init, d_lambda_init, d_lambda_min, ...
                                 it_newt_max, it_damp_max, tol, r_min, K_elast, Q, f, ...
                                 constitutive_matrix_builder, linear_system_solver)
%--------------------------------------------------------------------------
% init_phase_SSR_indirect_continuation performs the initial phase for the 
% indirect continuation methods. The purpose is to compute the 
% remaining input data needed by the continuation algorithms.
%
% In particular, given an initial strength reduction parameter lambda_init  
% and an initial increment d_lambda_init, this function computes:
%
%   lambda1 - The first value of lambda for which the Newton method converges.
%   lambda2 - The second value of lambda for which the Newton method converges.
%   U1      - The displacement field corresponding to lambda1.
%   U2      - The displacement field corresponding to the next parameter value lambda2.
%   omega1  - The value of the control parameter omega corresponding to lambda1.
%   omega2  - The value of omega corresponding to lambda2.
%   all_newton_its - An array containing the number of Newton iterations for each call.
%
% The function first solves the nonlinear system for lambda1, then 
% increases lambda by d_lambda and attempts to compute a new solution.
% If the solver fails for the second step, the increment d_lambda is halved 
% until a convergent solution is found. While omega2 is almost the same as 
% omega1, we replace U1, lambda1, omega1 with U2, lambda2, omega2, 
% respectively, set lambda2=lambda2+d_lambda1, and find new U2 and omega2.
%
% INPUTS:
%   lambda_init                 - Initial value of the strength reduction parameter.
%   d_lambda_init               - Initial increment in lambda.
%   d_lambda_min            - Minimal allowed increment in lambda.
%   it_newt_max             - Maximum number of Newton iterations.
%   it_damp_max             - Maximum number of damping iterations in Newton's method.
%   tol                     - Relative tolerance for Newton's solvers.
%   r_min                   - Basic regularization parameter for the stiffness matrix.
%   K_elast                 - Elastic stiffness matrix.
%   Q                       - Logical array restricting nodes with Dirichlet boundary conditions.
%   f                       - External load vector.
%   constitutive_matrix_builder - Object to build the constitutive force vector and tangent operator.
%   linear_system_solver    - Linear solver object (with deflation capabilities).
%
% OUTPUTS:
%   U1            - Displacement field corresponding to lambda1.
%   U2            - Displacement field corresponding to lambda2.
%   omega1        - Control parameter omega corresponding to lambda1.
%   omega2        - Control parameter omega corresponding to lambda2.
%   lambda1       - The updated value of lambda_init for which the Newton method converges.
%   lambda2       - The updated value of lambda for which the Newton method converges.
%   all_newton_its- A vector containing the number of Newton iterations from each call.
%
%--------------------------------------------------------------------------
%
% Initialization
%
n_n = size(Q, 2);                % Number of nodes.
dim = size(Q, 1);                % Spatial dimension (e.g., 2 or 3).
U_ini = zeros(dim, n_n);         % Initial displacement field (zero vector).

%
% First step of the continuation method.
%
fprintf('\n'); 
fprintf(' Step = %d  ', 1); 
fprintf('\n'); 

% Apply Davis' modifications to the material parameters using initial lambda.
constitutive_matrix_builder.reduction(lambda_init);
all_newton_its = [];         
% Newton's solver for initial lambda.
[U_it, flag_N, it_newt] = NEWTON.newton(U_ini, tol, it_newt_max, it_damp_max, r_min, ...
                       K_elast, Q, f, constitutive_matrix_builder, linear_system_solver);
all_newton_its = [all_newton_its it_newt];
% Test for convergence:
if flag_N == 1   % The solver was not successful.                      
     error('Initial choice of lambda seems to be too large.')
else           % The solver was successful.
     U1 = U_it;
     omega1 = f(Q)' * U1(Q);     
end

disp(['   lambda = ', num2str(lambda_init), ...
         ', d_lambda = ', num2str(d_lambda_init), ...
         ', omega = ', num2str(omega1)]); 

%    
% Second step of the continuation method.
%
fprintf('\n'); 
fprintf(' Step = %d  ', 2); 
fprintf('\n'); 

d_lambda = d_lambda_init;       % Set initial increment of lambda.
lambda1 = lambda_init;          % Set initial value of lambda
linear_system_solver.expand_deflation_basis(U_it(Q));
while true                   
                   
    % Update the parameter lambda.
    lambda_it = lambda1 + d_lambda;
                
    % Update strength parameters using lambda_it and Davis' approach.
    constitutive_matrix_builder.reduction(lambda_it);
               
    % Newton's solver.
    [U_it, flag_N, it_newt] = NEWTON.newton(U1, tol, it_newt_max, it_damp_max, r_min, ...
                       K_elast, Q, f, constitutive_matrix_builder, linear_system_solver);
    all_newton_its = [all_newton_its it_newt];
    % Test for convergence:
    if flag_N == 1  % The solver was not successful.
       % Decrease d_lambda.
       d_lambda = d_lambda / 2;
    else          % The solver was successful.
       U2 = U_it;
       omega2 = f(Q)' * U2(Q);     
       lambda2 = lambda_it;  
       if (abs(omega2-omega1)/omega1)<1e-5
           U1 = U2;
           lambda1 = lambda2;
           omega1 = omega2;
       else
           break 
       end
    end

    % Stopping criteria: if d_lambda becomes too small.
    if d_lambda < d_lambda_min
          error('It seems that the FoS is equal to lambda_init.')
    end
    if lambda1 > 10
          error('It seems that the FoS is greater than 10.')
    end
end     

disp(['   lambda = ', num2str(lambda2), ...
         ', d_lambda = ', num2str(lambda2 - lambda1), ...
         ', omega = ', num2str(omega1),...
         ', d_omega = ', num2str(omega2 - omega1)]);

end % function
