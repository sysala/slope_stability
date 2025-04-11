function [shur_prec_handle] = diag_prec_ICHOL(A, Q)
% diag_prec_ICHOL constructs a block-diagonal Schur complement preconditioner
% using incomplete Cholesky factorization.
%
% This function partitions the input matrix A into blocks based on the rows
% of the logical matrix Q. Each block is preconditioned using an incomplete
% Cholesky factorization (ichol) with a prescribed drop tolerance, and a
% solver function is defined for each block. The resulting function handle,
% shur_prec_handle, applies the block-diagonal preconditioner by solving the
% independent systems for each block and reassembling the solution.
%
% INPUTS:
%   A - Input matrix (typically the Schur complement of a larger system).
%       Its dimensions should match the number of true entries in Q (i.e.,
%       A is of size nnz(Q) x nnz(Q)).
%   Q - Logical matrix used to partition A into blocks. For example, if Q is a
%       d-by-N matrix with exactly one true entry per column, then each row of Q
%       defines a block.
%
% OUTPUT:
%   shur_prec_handle - Function handle that applies the block-diagonal
%                      preconditioner to a vector.
%
% NOTES:
%   - A small regularization (1e-2 times the diagonal) is added to each block
%     before performing the incomplete Cholesky factorization.
%   - The incomplete Cholesky parameters (e.g., drop tolerance) are set in the
%     'params' structure and can be adjusted as needed.
%

% Determine the number of blocks from the number of rows in Q.
dim = size(Q, 1);

% Preallocate cell arrays to store index vectors, diagonal blocks of A,
% and corresponding solver functions.
Q_cell = cell(dim, 1);
A_cell = cell(dim, 1);
L_cell = cell(dim, 1);
A_solver_cell = cell(dim, 1);

% Set incomplete Cholesky parameters.
params = struct('type', 'ict', 'droptol', 1e-3);

for i = 1:dim
    % Create a logical mask that retains only the i-th row of Q.
    temp = false(size(Q));
    temp(i, :) = Q(i, :);
    
    % Convert the mask into a vector of indices corresponding to the true entries.
    Q_cell{i} = temp(Q);
    
    % Extract the diagonal block of A corresponding to the i-th component.
    A_cell{i} = A(Q_cell{i}, Q_cell{i});
    
    % Regularize the block slightly and compute its incomplete Cholesky factor.
    L_cell{i} = ichol(A_cell{i} + 2e-2 * diag(diag(A_cell{i})), params);
    
    % Define the solver function for the current block.
    A_solver_cell{i} = @(x) L_cell{i}' \ (L_cell{i} \ x);
end

% Construct the overall block-diagonal preconditioner function handle.
shur_prec_handle = @(x) prec_apply_multi(x, Q_cell, A_solver_cell);
end

function res = prec_apply_multi(x, Q_cell, A_solver_cell)
% prec_apply_multi applies the block-diagonal incomplete Cholesky preconditioner.
%
% This function processes the input vector x by extracting subvectors
% corresponding to the index vectors in Q_cell, applying the corresponding
% solver (based on the incomplete Cholesky factorization) for each block,
% and reassembling the solution into a full vector.
%
% INPUTS:
%   x             - Input vector.
%   Q_cell        - Cell array of index vectors for each block (derived from Q).
%   A_solver_cell - Cell array of solver functions for each block.
%
% OUTPUT:
%   res - Output vector after applying the block-diagonal preconditioner.
%

% Initialize the result vector.
res = zeros(size(x));

% Number of blocks.
dim = numel(Q_cell);

for i = 1:dim
    % Extract the subvector corresponding to block i.
    b_tmp = x(Q_cell{i});
    
    % Apply the solver for block i.
    y_tmp = A_solver_cell{i}(b_tmp);
    
    % Place the computed solution into the corresponding indices.
    res(Q_cell{i}) = y_tmp;
end
end
