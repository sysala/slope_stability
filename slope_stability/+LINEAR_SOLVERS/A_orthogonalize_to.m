function [W_orth, norms] = A_orthogonalize_to(W, A, V, eps_add)
%--------------------------------------------------------------------------
% A_orthogonalize_to orthogonalizes the columns of V with respect to the A-inner
% product and an existing basis W.
%
% This function computes an orthogonal basis W_orth from the columns of V such that
% each new vector is orthogonal (in the A-inner product sense) to the columns of W
% and to those already computed in W_orth. The A-inner product of two vectors u and v 
% is defined as u'*(A*v). Vectors with an A-norm below eps_add are discarded.
%
% INPUTS:
%   W       - (Optional) Existing basis (matrix of size m x k) to be included 
%             in the orthogonalization process. If empty, only V is used.
%   A       - Symmetric positive definite matrix defining the A-inner product.
%   V       - Matrix (size m x n) whose columns are to be orthogonalized.
%   eps_add - Tolerance for the A-norm of a vector to be considered nonzero.
%
% OUTPUTS:
%   W_orth  - Matrix whose columns form an A-orthogonal basis extracted from V.
%   norms   - Vector of scaling factors for the corresponding columns of W_orth.
%             (A positive value indicates a standard normalization; a negative value
%             indicates that the computed A-norm was negative before taking the absolute.)
%
%--------------------------------------------------------------------------
[m, n] = size(V);

W_orth = zeros(m, n);
norms = ones(n, 1);
mask_add = true(1, n);

for i = 1:n
    temp = V(:, i);
    A_temp = A * temp;
    % Remove components in the direction of the existing basis W (if provided)
    if ~isempty(W)
        temp = temp - W * ((A_temp)' * W)';
    end
    % Remove components in the direction of already computed orthogonal vectors.
    if i > 1
        temp = temp - (W_orth(:, 1:(i-1)) .* norms(1:(i-1))') * ((A_temp)' * W_orth(:, 1:(i-1)))';
    end

    norm_temp = temp' * A_temp;
    if abs(norm_temp) > eps_add
        if norm_temp > 0
            norm_temp = sqrt(norm_temp);
        else
            norm_temp = sqrt(abs(norm_temp));
            norms(i) = -1;
        end
        W_orth(:, i) = temp / norm_temp;
    else
        mask_add(i) = 0;
    end
end

W_orth = W_orth(:, mask_add);
norms = norms(mask_add);
end
