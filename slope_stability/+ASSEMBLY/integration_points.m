function [transformed_points] = integration_points(elem,coord,Xi)

% Assuming the following variables are defined:
% elem  : 10 x m matrix, each column contains node indices for an element
% coord : 3 x k matrix, each column contains the (x, y, z) coordinates of a node
% Xi    : 3 x n matrix, each column contains coordinates in the reference tetrahedron

% Number of elements (m) and number of points (n)
m = size(elem, 2);
n = size(Xi, 2);

% Extract the first four nodes of each element for affine mapping
% v0, v1, v2, v3 correspond to the vertices of the tetrahedron
v0 = coord(:, elem(1, :)); % 3 x m
v1 = coord(:, elem(2, :)); % 3 x m
v2 = coord(:, elem(3, :)); % 3 x m
v3 = coord(:, elem(4, :)); % 3 x m

% Compute the edge vectors relative to v0
d1 = v1 - v0; % 3 x m
d2 = v2 - v0; % 3 x m
d3 = v3 - v0; % 3 x m

% Reshape for broadcasting to perform vectorized operations
% Adding a singleton dimension for points
v0_exp  = reshape(v0, 3, m, 1); % 3 x m x 1
d1_exp  = reshape(d1, 3, m, 1); % 3 x m x 1
d2_exp  = reshape(d2, 3, m, 1); % 3 x m x 1
d3_exp  = reshape(d3, 3, m, 1); % 3 x m x 1

% Reshape Xi components to align with elements
xi    = reshape(Xi(1, :), 1, 1, n); % 1 x 1 x n
eta   = reshape(Xi(2, :), 1, 1, n); % 1 x 1 x n
zeta  = reshape(Xi(3, :), 1, 1, n); % 1 x 1 x n

% Perform the affine transformation
% x_all will be a 3 x m x n array containing all transformed points
x_all = v0_exp + d1_exp .* xi + d2_exp .* eta + d3_exp .* zeta; % 3 x m x n

% Rearrange the dimensions to have elements first, then points
x_all = permute(x_all, [1, 3, 2]); % 3 x n x m

% Reshape to a 3 x (n*m) matrix with the desired ordering
transformed_points = reshape(x_all, 3, m * n); % 3 x (m*n)