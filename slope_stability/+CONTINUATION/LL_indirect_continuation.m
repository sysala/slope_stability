function [U, t_hist, omega_hist, U_max_hist] = LL_indirect_continuation( ...
    d_omega_ini, d_t_min, step_max, omega_max, ...
    it_newt_max, it_damp_max, tol, r_min, K_elast, U_elast, Q, f, ...
    constitutive_matrix_builder, linear_system_solver)
%--------------------------------------------------------------------------
% LL_indirect_continuation solves a limit load problem by an indirect 
% continuation technique. The control parameter "omega" (representing the
% work of external forces) is increased adaptively, and for each omega we 
% solve the nonlinear system:
%
%       F(U) = t * f,    and    f' * U = omega
%
% where U is the displacement vector, and t is the load parameter whose
% limit value t* is unknown. As omega → +∞, t → t*.
%
% The increment d_omega may be increased/decreased adaptively based on 
% convergence success/failure and the observed increments in t.
%
% INPUT ARGUMENTS:
%   d_omega_ini - Initial increment of omega.
%   d_t_min     - Critical lower bound for the increment d_t of the load factor t.
%   step_max    - Maximum number of load steps.
%   omega_max   - Maximum value of the control parameter omega.
%   it_newt_max - Maximum number of Newton iterations per step.
%   it_damp_max - Maximum number of damping steps for each Newton iteration.
%   tol         - Tolerance for stopping Newton's method.
%   r_min       - Basic regularization parameter for the stiffness matrix.
%   K_elast     - Elastic stiffness matrix.
%   U_elast     - Elastic displacement field (initial guess).
%   Q           - Logical array (3 x n_n) restricting Dirichlet boundary DOFs.
%   f           - External load vector.
%   constitutive_matrix_builder - Object to build the constitutive force vector 
%                                 and tangent operator.
%   linear_system_solver - Linear solver object (with deflation capability).
%
% OUTPUT:
%   U          - Final displacement field at the last converged omega.
%   t_hist     - History of the load parameter t for each successful step.
%   omega_hist - History of the control parameter omega for each successful step.
%   U_max_hist - History of the maximum displacement per step.
%
%--------------------------------------------------------------------------

% Pre-allocate history arrays for speed.
omega_hist  = zeros(1, step_max);
t_hist      = zeros(1, step_max);
U_max_hist  = zeros(1, step_max);

% Initialization of continuation parameters.
omega       = 0;              % Current control parameter.
omega_old   = 0;              % Previous value of omega.
d_omega     = d_omega_ini;    % Current increment of omega.
t           = 0;              % Current load parameter.
U           = U_elast;        % Current displacement (start from elastic solution).
U_old       = U;              % Previous displacement.
U_ini       = U;              % Initial guess for Newton in each step.

% For limiting repeated halving of d_omega.
n_omega     = 0;              % Number of consecutive reductions in d_omega.
n_omega_max = 5;              % Maximum number of reductions allowed.

% Start counting the continuation steps.
step = 1;

% Main continuation loop.
while true
    
    fprintf('\nStep = %d\n', step+1);
    
    % Try to increase the control parameter by d_omega.
    omega_it = omega + d_omega;
    if omega_it > omega_max
        % Do not exceed the specified maximum.
        omega_it = omega_max;
        d_omega  = omega_it - omega; 
    end

    % Simple initial guess for this step.
    if step > 1
        U_ini = U; 
    end

    %----------------------------------------------------------------------
    %   Call Newton's method (with damping, regularization, etc.)
    %   to attempt solving for U_it and t_it given the new omega_it.
    %   "flag == 0" indicates success; "flag == 1" indicates failure.
    %----------------------------------------------------------------------
    [U_it, t_it, flag, newton_its] = NEWTON.newton_ind_LL( ...
        U_ini, t, omega_it, ...
        it_newt_max, it_damp_max, tol, r_min, K_elast, Q, f, ...
        constitutive_matrix_builder, linear_system_solver.copy());
    
    % Check convergence of the Newton solver.
    if flag == 1  % Solver was not successful.
        d_omega = d_omega / 2;
        n_omega = n_omega + 1;
    else
        % Solver was successful: update solution and continuation parameters.
        step       = step + 1;
        U_old      = U;
        U          = U_it;
        linear_system_solver.expand_deflation_basis(U(Q));
        omega_old  = omega;
        omega      = omega_it;
        
        d_t        = t_it - t;
        t          = t_it;
        n_omega    = 0;  % Reset consecutive reduction counter.

        % Store history.
        omega_hist(step) = omega;
        t_hist(step)     = t;
        
        % Compute maximum nodal displacement for monitoring.
        U_total = sqrt(sum(U.^2));
        U_max   = max(U_total);
        U_max_hist(step) = U_max;
        
        % Display iteration information.
        fprintf('   t = %g, d_t = %g, omega = %g, d_omega = %g, U_max = %g\n', ...
                t, d_t, omega, d_omega, U_max);
        
        % If the scaled increase in lambda is minimal, stop continuation.
        if (d_t < d_t_min)
            disp('Minimal increment of lambda was achieved.');
            break;
        end
        
        % Criterion for the change of d_omega.
        if (newton_its < 20) && (step > 2) && ...
           ((t_hist(step) - t_hist(step-1)) < 0.9 * (t_hist(step-1) - t_hist(step-2)))
            d_omega = 2 * d_omega;
        end
    end
    
    % Stopping criteria.
    if n_omega >= n_omega_max
        disp('It is impossible to increase omega further (max reductions reached).');
        break;
    end
    if omega >= omega_max
        disp('Maximal omega was achieved.');
        break;
    end
    % if U_max > U_max_stop
    %     disp('Maximal displacement achieved.');
    %     break;
    % end
    if step >= step_max
        disp('Maximal number of steps was reached.');
        break;
    end

end % while

% Clip the output arrays to the actual number of successful steps.
t_hist     = t_hist(1:step);
omega_hist = omega_hist(1:step);
U_max_hist = U_max_hist(1:step);

end % function
