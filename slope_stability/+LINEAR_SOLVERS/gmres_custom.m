function [x, flag, relres, iter, resvec] = gmres_custom(A, b, restart, tol, maxit, x0)
% GMRES_CUSTOM Solve Ax=b using GMRES algorithm with restart
%
% Inputs:
%   A       - Matrix or function handle for matrix-vector multiplication
%   b       - Right-hand side vector
%   restart - Number of iterations before restart
%   tol     - Tolerance for convergence
%   maxit   - Maximum number of iterations (outer iterations)
%   x0      - Initial guess
%
% Outputs:
%   x       - Solution vector
%   flag    - Convergence flag (0 if converged, 1 otherwise)
%   relres  - Relative residual
%   iter    - Number of iterations performed
%   resvec  - Vector containing residual norms

if nargin < 6
    x0 = zeros(size(b));
end

n = length(b);
x = x0;
r = b - A*x;
beta = norm(r);
resvec = beta;
relres = beta / norm(b);

if relres < tol
    flag = 0;
    iter = 0;
    return;
end

flag = 1;
iter = 0;

for outer = 1:maxit
    V = zeros(n, restart+1);
    H = zeros(restart+1, restart);
    cs = zeros(restart, 1);
    sn = zeros(restart, 1);

    V(:,1) = r / beta;
    s = beta * eye(restart+1, 1);

    for j = 1:restart
        iter = iter + 1;

        % Arnoldi process
        w = A * V(:,j);
        for i = 1:j
            H(i,j) = V(:,i)' * w;
            w = w - H(i,j)*V(:,i);
        end
        H(j+1,j) = norm(w);
        if H(j+1,j) ~= 0
            V(:,j+1) = w / H(j+1,j);
        end

        % Apply Givens rotations
        for i = 1:j-1
            temp = cs(i)*H(i,j) + sn(i)*H(i+1,j);
            H(i+1,j) = -sn(i)*H(i,j) + cs(i)*H(i+1,j);
            H(i,j) = temp;
        end

        % Compute and apply new Givens rotation
        [cs(j), sn(j)] = givens_rotation(H(j,j), H(j+1,j));
        temp = cs(j)*s(j);
        s(j+1) = -sn(j)*s(j);
        s(j) = temp;

        % Eliminate H(j+1,j)
        H(j,j) = cs(j)*H(j,j) + sn(j)*H(j+1,j);
        H(j+1,j) = 0;

        % Update residual
        res = abs(s(j+1));
        resvec = [resvec; res];
        relres = res / norm(b);

        if relres < tol
            y = H(1:j,1:j) \ s(1:j);
            x = x + V(:,1:j)*y;
            flag = 0;
            return;
        end
    end

    % Update solution and residual
    y = H(1:restart, 1:restart) \ s(1:restart);
    x = x + V(:,1:restart)*y;
    r = b - A*x;
    beta = norm(r);
    relres = beta / norm(b);

    if relres < tol
        flag = 0;
        return;
    end

    resvec = [resvec; beta];
end

end

function [c,s] = givens_rotation(a,b)
% Compute Givens rotation parameters c and s
if b == 0
    c = 1;
    s = 0;
else
    if abs(b) > abs(a)
        temp = a / b;
        s = 1 / sqrt(1 + temp^2);
        c = temp * s;
    else
        temp = b / a;
        c = 1 / sqrt(1 + temp^2);
        s = temp * c;
    end
end
end
