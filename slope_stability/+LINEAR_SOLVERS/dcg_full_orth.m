function [x, iter, resvec] = dcg_full_orth(A, b, M, W, maxits, tol, x0)
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

r=b-A*x;
b_norm=norm(b);
res=norm(r)/b_norm;
if res<tol || maxits==0 || b_norm==0
    resvec=res;
    iter=0;
    return
end
z=M(r);
p=P(z);

gamma_old=dot(r,z);
resvec=zeros(maxits+1,1);
resvec(1)=res;

all_p = cell(maxits,1);
all_Ap = cell(maxits,1);

for j=1:maxits
    s=A*p;
    p_A_norm = sqrt(dot(s, p));
    all_p{j} = p/p_A_norm;
    all_Ap{j} = s/p_A_norm;
    
    alpha=gamma_old/dot(s,p);
    x=x+alpha*p;
    r=r-alpha*s;
    res = norm(r)/b_norm;
    resvec(j+1)=res;
    if res<tol
        break;
    end
    z=M(r);
    gamma_new=dot(r,z);
    for k=1:j
        dot_tmp = dot(z,all_Ap{k});
        z = z - dot_tmp*all_p{k};
    end
    % beta=gamma_new/gamma_old;
    % p=P(z)+beta*p;
    p=P(z);
    gamma_old=gamma_new;
end

resvec=resvec(1:j+1);
iter=j;
end
