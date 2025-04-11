function [W_orth, norms] = A_orthogonalize(W, A, eps_add)
%--------------------------------------------------------------------------
% A_orthogonalize computes an A-orthogonal basis from the columns of W.
%
% This function orthogonalizes the columns of W with respect to the A-inner
% product. The A-inner product between two vectors u and v is defined as
% u'*(A*v). The columns of W are processed in reverse order, and only those
% with an A-norm above the tolerance eps_add are retained.
%
% INPUTS:
%   W      - Matrix whose columns are to be A-orthogonalized (size m x n).
%   A      - Symmetric positive definite matrix defining the A-inner product.
%   eps_add- Tolerance for the A-norm; vectors with norm below eps_add are discarded.
%
% OUTPUTS:
%   W_orth - Matrix whose columns form an A-orthogonal basis (subset of columns of W).
%   norms  - Vector containing scaling factors for the corresponding columns in W_orth.
%            A positive value indicates standard normalization; a negative value
%            indicates that the computed A-norm was originally negative before taking
%            the absolute value.
%
%--------------------------------------------------------------------------
[m, n] = size(W);

W_orth = zeros(m, n);
norms = ones(n, 1);
mask_add = false(1, n);

% Process columns in reverse order.
for i = 1:n
    % Get the i-th vector from the end of W.
    temp = W(:, n - i + 1);
    A_temp = A * temp;
    % Remove the contributions from previously computed orthogonal vectors.
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
        mask_add(i) = true;
    end
end

% Retain only the columns that were accepted.
W_orth = W_orth(:, mask_add);
norms = norms(mask_add);
% Reverse the order so that the columns correspond to the original order.
W_orth = W_orth(:, end:-1:1);
norms = norms(end:-1:1);
end
