function [x, nit, err, resvec] = flexible_GMRES_deflate(A, b, m, tol, maxit, Prec_fct, W)
%--------------------------------------------------------------------------
% flexible_GMRES_deflate solves the linear system Ax = b using Flexible GMRES
% with deflation and preconditioning.
%
% This function implements a Flexible GMRES algorithm with a restart parameter
% m and deflation. If a deflation basis W is provided, the method constructs 
% coarse space correction operators. Preconditioning is applied via the function 
% Prec_fct. The algorithm uses Givens rotations for the least-squares problem.
%
% INPUTS:
%   A       - Matrix A (numeric or function handle).
%   b       - Right-hand side vector.
%   m       - Dimension of the Krylov subspace (restart parameter).
%   tol     - Tolerance for the relative residual.
%   maxit   - Maximum number of iterations.
%   Prec_fct- Preconditioner function (matrix or function handle) to be applied.
%   W       - Deflation basis matrix (columns are deflation vectors). If empty,
%             no deflation is performed.
%
% OUTPUTS:
%   x       - Computed solution vector.
%   nit     - Number of iterations performed.
%   err     - Final relative residual.
%   resvec  - Vector of residual norms at each iteration.
%
% Reference:
%   Saad, "A flexible inner-outer preconditioned GMRES algorithm", SIAM 1993
%--------------------------------------------------------------------------
    
% If no deflation basis is provided, define identity operators.
if isempty(W)
    Proj_fct = @(x) x;
    Q_coarse_solve = @(x) 0 * x;
else
    % Define coarse solve and projection operators based on deflation basis W.
    Q_coarse_solve = @(x) W * (W' * x);
    Proj_fct = @(x) x - W * (W' * (A * x));
    %Q_coarse_solve = @(x) 0 * x;
    %Proj_fct = @(x) x;
end

% Get problem size and initialize residual storage.
N = size(b, 1);
resvec = ones(1, maxit);

% Wrap A in a function handle if necessary.
Afct = @(u) A * u;

% Compute an initial guess using the coarse solve.
x0 = Q_coarse_solve(b);

% Initialize variables for FGMRES.
normb = norm(b);
cs = zeros(m, 1);
sn = zeros(m, 1);
e1 = zeros(N, 1);
e1(1) = 1.0;
nit = 1;

% Begin the FGMRES iteration loop.
while nit < maxit
    % Initialize the Hessenberg matrix and the initial residual.
    H = zeros(m+1, m);
    r0 = b - Afct(x0);
    beta = norm(r0);
    V = zeros(N, m+1);
    V(:, 1) = (1 / beta) * r0;
    s = beta * e1;
    Z = zeros(N, m);
    
    % Arnoldi process: construct the Krylov subspace.
    for i = 1:m
        % Apply the preconditioner and deflation.
        tmp = Prec_fct(V(:, i));
        Z(:, i) = Proj_fct(tmp);
        w = Afct(Z(:, i));
        
        % Modified Gram-Schmidt orthogonalization.
        for k = 1:i
            H(k, i) = w' * V(:, k);
            w = w - H(k, i) * V(:, k);
        end
        
        H(i+1, i) = norm(w);
        V(:, i+1) = (1 / H(i+1, i)) * w;
        
        % Apply previously computed Givens rotations.
        for k = 1:i-1
            temp = cs(k) * H(k, i) + sn(k) * H(k+1, i);
            H(k+1, i) = -sn(k) * H(k, i) + cs(k) * H(k+1, i);
            H(k, i) = temp;
        end
        
        % Compute and apply a new Givens rotation.
        [cs(i), sn(i)] = givens_rotmat(H(i, i), H(i+1, i));
        temp = cs(i) * s(i);
        s(i+1) = -sn(i) * s(i);
        s(i) = temp;
        H(i, i) = cs(i) * H(i, i) + sn(i) * H(i+1, i);
        H(i+1, i) = 0.0;
        err = abs(s(i+1)) / normb;
        
        resvec(nit) = err;  % Store current residual norm.
        
        % Check for convergence.
        if err <= tol
            y = H(1:i, 1:i) \ s(1:i);
            x = x0 + Z(:, 1:i) * y;
            resvec = resvec(1:nit);
            return
        end
        
        if nit == maxit
            break;
        end
        
        nit = nit + 1;
    end
    
    % Update the solution after m Arnoldi iterations.
    if err <= tol
        break;
    end
    y = H(1:m, 1:m) \ s(1:m);
    x = x0 + Z(:, 1:m) * y;
    r = b - Afct(x);
    s(i+1) = norm(r);
    err = s(i+1) / normb;
    
    if err <= tol
        break;
    end
    
    % Update initial guess for next cycle.
    x0 = x;
    resvec = resvec(1:min(end, nit));
end

%--------------------------------------------------------------------------
% Helper function: givens_rotmat
% Computes the Givens rotation parameters (c,s) to eliminate b.
function [c, s] = givens_rotmat(a, b)
    if b == 0.0
        c = 1.0;
        s = 0.0;
    elseif abs(b) > abs(a)
        temp = a / b;
        s = 1.0 / sqrt(1.0 + temp^2);
        c = temp * s;
    else
        temp = b / a;
        c = 1.0 / sqrt(1.0 + temp^2);
        s = temp * c;
    end
end

end
