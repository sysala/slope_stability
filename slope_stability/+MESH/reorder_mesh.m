function [coord, elem, surf, Q, perm, iperm] = reorder_mesh(coord, elem, surf, Q)
% MESH.reorder_mesh  Reorder mesh nodes to reduce matrix bandwidth.
%
%   [coord, elem, surf, Q, perm, iperm] = MESH.reorder_mesh(coord, elem, surf, Q)
%
%   Uses Octave's built-in Reverse Cuthill-McKee (symrcm) ordering to
%   reduce the bandwidth of the stiffness matrix, improving cache locality
%   for iterative solvers and preconditioner applications.
%
%   Inputs:
%     coord   – (dim, n_n)  node coordinates (dim = 2 or 3)
%     elem    – (n_loc, n_e)  element connectivity (1-based)
%     surf    – (n_sf, n_s) surface/edge connectivity (1-based), or []
%     Q       – (dim, n_n)  logical free-DOF mask
%
%   Outputs:
%     coord, elem, surf, Q  – reordered mesh data
%     perm   – (n_n, 1)  old-to-new mapping:  perm(old) = new
%     iperm  – (n_n, 1)  new-to-old mapping:  iperm(new) = old
%
%   The returned perm/iperm can be used to reorder other node-indexed arrays:
%     Q_w_new  = Q_w(iperm);        % or equivalently: Q_w_new(perm) = Q_w;
%     pw_D_new = pw_D(iperm);
%
%   After reordering, all assembly functions (elastic_stiffness_matrix_*,
%   vector_volume_*, seepage, etc.) should be called with the reordered data.

n_n   = size(coord, 2);
n_e   = size(elem, 2);
n_loc = size(elem, 1);

% ---- Build node adjacency graph from element connectivity ----
% Build sparse incidence E (n_n x n_e), then adjacency A = E*E' - I.
elem_d = double(elem);
row_idx = elem_d(:);
col_idx = repmat(1:n_e, n_loc, 1);  col_idx = col_idx(:);
E = sparse(row_idx, col_idx, 1, n_n, n_e);
A = spones(E * E');
A = A - spdiags(diag(A), 0, n_n, n_n);

% ---- Reverse Cuthill-McKee ordering (built-in) ----
iperm = symrcm(A)';    % new-to-old
perm  = zeros(n_n, 1); % old-to-new
perm(iperm) = (1:n_n)';

% ---- Apply permutation ----
coord = coord(:, iperm);
elem  = int32(perm(elem));
if ~isempty(surf)
    surf = int32(perm(surf));
end
Q = Q(:, iperm);

% Report bandwidth reduction
bw_old = max(abs(diff(find(A))));
A_new  = A(iperm, iperm);
bw_new = max(abs(diff(find(A_new))));
fprintf('  Mesh reordering (symrcm): bandwidth %d -> %d (%.1fx reduction)\n', ...
    bw_old, bw_new, bw_old / max(bw_new, 1));
end
