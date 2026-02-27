function [Z_free, labels, Z_full] = near_null_space_elasticity_3D(coord, Q, center_coordinates)
%--------------------------------------------------------------------------
% near_null_space_elasticity_3D builds near-null-space modes for 3D linear
% elasticity from mesh coordinates and displacement boundary mask.
%
% The function constructs the 6 rigid-body modes:
%   3 translations + 3 rotations, then restricts them to free DOFs (Q==true).
% The output basis is Euclidean-orthonormalized with a stable Gram-Schmidt
% process and near-zero dependent modes are automatically removed.
%
% INPUTS:
%   coord  - Nodal coordinates, size (3, n_n)
%   Q      - Logical free-DOF mask, size (3, n_n)
%   center_coordinates (optional, default=true)
%          - If true, rotations use coordinates shifted by centroid.
%
% OUTPUTS:
%   Z_free - Near-null-space basis on free DOFs, size (nnz(Q), k), k<=6
%   labels - Cell array with names of retained modes (length k)
%   Z_full - Same basis in full DOFs, size (3*n_n, k)
%--------------------------------------------------------------------------

if nargin < 3 || isempty(center_coordinates)
    center_coordinates = true;
end

if size(coord, 1) ~= 3
    error('coord must have size (3, n_n).');
end
if ~isequal(size(Q), size(coord))
    error('Q must have the same size as coord.');
end

n_n = size(coord, 2);
free_mask = logical(Q(:));

x = coord(1, :);
y = coord(2, :);
z = coord(3, :);
if center_coordinates
    x = x - mean(x);
    y = y - mean(y);
    z = z - mean(z);
end

mode_names = {'tx', 'ty', 'tz', 'rx', 'ry', 'rz'};
modes_full = zeros(3 * n_n, 6);

% Translation in x, y, z.
U = zeros(3, n_n);
U(1, :) = 1;
modes_full(:, 1) = U(:);
U = zeros(3, n_n);
U(2, :) = 1;
modes_full(:, 2) = U(:);
U = zeros(3, n_n);
U(3, :) = 1;
modes_full(:, 3) = U(:);

% Rotation modes: u = omega x r.
U = zeros(3, n_n);        % around x-axis
U(2, :) = -z;
U(3, :) =  y;
modes_full(:, 4) = U(:);

U = zeros(3, n_n);        % around y-axis
U(1, :) =  z;
U(3, :) = -x;
modes_full(:, 5) = U(:);

U = zeros(3, n_n);        % around z-axis
U(1, :) = -y;
U(2, :) =  x;
modes_full(:, 6) = U(:);

% Restrict to free DOFs.
Z_candidates = modes_full(free_mask, :);

% Stable orthonormalization with dependency filtering (preserve mode order).
tol = 1e-12 * sqrt(max(1, size(Z_candidates, 1)));
Z_free = zeros(size(Z_candidates, 1), 6);
Z_full = zeros(3 * n_n, 6);
labels = {};
n_kept = 0;

for i = 1:6
    v_free = Z_candidates(:, i);
    v_full = modes_full(:, i);

    if n_kept > 0
        coeff = Z_free(:, 1:n_kept)' * v_free;
        v_free = v_free - Z_free(:, 1:n_kept) * coeff;
        v_full = v_full - Z_full(:, 1:n_kept) * coeff;
    end

    nv = norm(v_free);
    if nv > tol
        n_kept = n_kept + 1;
        Z_free(:, n_kept) = v_free / nv;
        Z_full(:, n_kept) = v_full / nv;
        labels{n_kept, 1} = mode_names{i}; %#ok<AGROW>
    end
end

Z_free = Z_free(:, 1:n_kept);
Z_full = Z_full(:, 1:n_kept);

end
