function [x, iter, resvec, mresvec, errvec] = dcg_M(A, b, M, W, maxits, tol, x0, x_exact)
%DEFLATEDCG  Deflated Conjugate Gradient solver.
%
%   [x, iters, res_hist] = DCG(A, b, M, W, maxits, tol, x0) solves the
%   symmetric positive definite linear system A*x = b using a deflated version of
%   the Conjugate Gradient method. When a deflation basis matrix W is provided (non-empty),
%   a coarse solve and projection operator are used to remove components associated with
%   the deflation subspace.
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
%       When deflation is activated (W is non-empty), the algorithm uses the coarse
%       solve operator
%           Q_coarse_solve = @(x) W * (W' * x),
%       and the projection operator
%           Proj_fct = @(x) x - W * (W' * (A * x)).
%       These operators default to the identity if W is empty.
%
%   See also CG, PCG.

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
    P = @(x) x;
else
    % Define coarse solve and projection operators based on deflation basis W.
    Q_coarse_solve = @(x) W * (W' * x);
    P = @(x) x - W * (W' * (A * x));
    x = Q_coarse_solve(b);
end
error_instead = 1;
if isempty(x_exact)
    error_instead = 0;
end
r=b-A*x;
b_norm=norm(b);
Mb_norm=norm(M(b));

resvec = zeros(maxits+1,1);
mresvec = zeros(maxits+1,1);
errvec = zeros(maxits+1,1);

err=0;
if error_instead
    err=norm(x-x_exact);
end

mres=norm(M(A*x-b))/Mb_norm;
res = norm(r)/b_norm;

resvec(1)=res;
mresvec(1)=mres;
errvec(1)=err;

if res<tol || maxits==0 || b_norm==0
    resvec=res;
    mresvec=mres;
    errvec=err;
    iter=0;
    return
end
z=M(r);
p=P(z);

gamma_old=dot(r,z);


for j=1:maxits
    s=A*p;
    alpha=gamma_old/dot(s,p);
    x=x+alpha*p;
    r=r-alpha*s;


    if error_instead
        err=norm(x-x_exact);
    end
    mres=norm(M(A*x-b))/Mb_norm;
    res = norm(A*x - b)/b_norm;
    resvec(j+1)=res;
    mresvec(j+1)=mres;
    errvec(j+1)=err;

    
    if res<tol
        break;
    end
    z=M(r);
    gamma_new=dot(r,z);
    beta=gamma_new/gamma_old;
    p=P(z)+beta*p;
    gamma_old=gamma_new;
end

resvec=resvec(1:j+1);
mresvec=mresvec(1:j+1);
errvec=errvec(1:j+1);
iter=j;
end
