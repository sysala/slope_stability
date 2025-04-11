function [shur_prec_handle] = diag_prec_AGMG(A, Q)
% diag_prec_AGMG constructs a block-diagonal Schur complement preconditioner
% using the AGMG solver.
%
% This function partitions the input matrix A into blocks based on the
% logical matrix Q. Each block corresponds to one row of Q. For each block,
% the AGMG multigrid solver is initialized with specific parameters and a
% solver function is defined to apply AGMG to a vector. The returned function
% handle, shur_prec_handle, applies the block-diagonal preconditioner by solving
% the independent systems for each block and assembling the results.
%
% INPUTS:
%   A - Input matrix (typically the Schur complement of a larger system).
%       Its dimensions should match the number of true entries in Q (i.e.,
%       A is of size nnz(Q) x nnz(Q)).
%   Q - Logical matrix used to partition A into blocks. For example, if Q is a
%       3-by-N matrix with exactly one true entry per column, then each row of Q
%       defines a block.
%
% OUTPUT:
%   shur_prec_handle - Function handle that applies the block-diagonal
%                      preconditioner to a vector.
%
% NOTES:
%   - AGMG must be available in your MATLAB path.
%   - AGMG initialization and solver parameters are hard-coded and may be
%     adjusted as necessary.
%

% Determine the number of blocks from the number of rows in Q.
dim = size(Q,1);

% Preallocate cell arrays to store index vectors, diagonal blocks of A,
% and the corresponding solver functions.
Q_cell = cell(dim,1);
A_cell = cell(dim,1);
A_solver_cell = cell(dim,1);

for i = 1:dim
    % Create a logical mask that retains only the i-th row of Q.
    temp = false(size(Q));
    temp(i,:) = Q(i,:);
    
    % Convert the logical mask into a vector of indices (using column-major
    % order) corresponding to the true entries.
    Q_cell{i} = temp(Q);
    
    % Extract the diagonal block from A corresponding to the current indices.
    A_cell{i} = A(Q_cell{i}, Q_cell{i});
    
    % Initialize AGMG for the current block with specific parameters.
    agmg(A_cell{i}, [], 1, [], [], -1, [], 1, 8, 4, i);
    
    % Define the solver function for the current block.
    A_solver_cell{i} = @(x) agmg(A_cell{i}, x, 1, 1e-6, 10, -1, [], 3, 8, 4, i);
end

% Construct the overall block-diagonal preconditioner as a function handle.
shur_prec_handle = @(x) prec_apply(x, Q_cell, A_solver_cell);
end

function [res] = prec_apply(x, Q_cell, A_solver_cell)
% prec_apply applies the block-diagonal AGMG preconditioner to a vector x.
%
% This function processes the input vector x by extracting subvectors
% corresponding to the indices in each block (stored in Q_cell), applying
% the respective AGMG solver (from A_solver_cell) to each subvector, and then
% reassembling the results into a full output vector.
%
% INPUTS:
%   x             - Input vector.
%   Q_cell        - Cell array of index vectors for each block (derived from Q).
%   A_solver_cell - Cell array of solver functions for each corresponding block.
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
    
    % Apply the AGMG solver for block i.
    y_tmp = A_solver_cell{i}(b_tmp);
    
    % Place the computed solution back into the appropriate locations.
    res(Q_cell{i}) = y_tmp;
end
end
