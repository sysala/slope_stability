function [x, iters, res_hist] = dfgmres_solver(A, b, M, W, maxits, tol, x0)
%DFGMRES  Deflated Flexible GMRES iterative solver.
%
%   [x, iters, res_hist] = DFGMRES(A, b, M, W, maxits, tol, x0) solves the
%   linear system A*x = b using a deflated version of the flexible GMRES method.
%
%   The algorithm incorporates deflation through a coarse space provided via the
%   deflation basis matrix W. When W is non-empty, an initial coarse solve is applied
%   along with a projection operator to accelerate convergence. If W is empty, the
%   method reverts to the standard FGMRES approach.
%
%   Inputs:
%       A      - Function handle for the matrix-vector product, i.e., A(x)
%       b      - Right-hand side vector
%       M      - Function handle for the preconditioner, i.e., M(x)
%       W      - Deflation basis matrix. If empty, no deflation is performed.
%       maxits - Maximum number of iterations (default: 1000)
%       tol    - Relative residual tolerance (default: 1e-6)
%       x0     - Initial guess (default: zero vector)
%
%   Outputs:
%       x       - Computed solution vector
%       iters   - Number of iterations performed
%       res_hist- History of the relative residual norms
%
%   Note:
%       The final solution is reconstructed as
%           x = x0 + sum_{j=1}^{iters} y(j) * w_j,
%       where y is the solution of the small least-squares problem arising from the Arnoldi process,
%       and w_j are the preconditioned (and projected) vectors.
%
%   See also GMRES, PCG.

if isempty(maxits)
    maxits = 1000;
end
if isempty(tol)
    tol = 1e-6;
end
if isempty(x0)
    x = zeros(size(b));
else
    x = x0;
end
if isempty(W)
    Proj_fct = @(x) x;
else
    % Define coarse solve and projection operators based on deflation basis W.
    Q_coarse_solve = @(x) W * (W' * x);
    Proj_fct = @(x) x - W * (W' * (A * x));
    x = Q_coarse_solve(b);
end


n = length(b);
r0 = b - A*x;
b_norm = norm(b);
if b_norm == 0, b_norm = 1; end
res_norm = norm(r0);

% Preallocate the residual history vector
res_hist = zeros(maxits+1, 1);
res_hist(1) = res_norm / b_norm;

if res_hist(1) < tol
    iters = 0;
    res_hist = res_hist(1);
    return;
end

% Preallocate arrays for Arnoldi vectors, Hessenberg matrix, and preconditioned vectors
V = zeros(n, maxits+1);
H = zeros(maxits+1, maxits);
W = zeros(n, maxits);

V(:,1) = r0 / res_norm;
g = zeros(maxits+1, 1);
g(1) = res_norm;

for j = 1:maxits
    % Apply preconditioner to the j-th Krylov basis vector
    w = M(V(:,j));
    w = Proj_fct(w);
    % Compute A*w
    u = A*w;
    
    % Arnoldi process
    for i = 1:j
        H(i,j) = V(:,i)' * u;
        u = u - H(i,j) * V(:,i);
    end
    H(j+1,j) = norm(u);
    W(:,j) = w;  % store the preconditioned vector

    % Check for happy breakdown
    if H(j+1,j) < 1e-14
        iters = j;
        break;
    end

    V(:,j+1) = u / H(j+1,j);

    % Solve the least-squares problem
    y = H(1:j+1, 1:j) \ g(1:j+1);

    % Compute the true residual norm explicitly
    res_norm = norm(g(1:j+1) - H(1:j+1, 1:j) * y);
    res_hist(j+1) = res_norm / b_norm;

    if res_hist(j+1) < tol
        iters = j;
        break;
    end
    iters = j;  % update iteration count
end

% Trim the residual history to the iterations performed
res_hist = res_hist(1:iters+1);

% Reconstruct the final solution: x = x0 + sum_{j=1}^{iters} y(j)*W(:,j)
x = x + W(:,1:iters) * y;

end
