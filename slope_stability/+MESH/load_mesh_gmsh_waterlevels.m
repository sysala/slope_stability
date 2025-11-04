function [coord, elem, surf, Q, material, triangle_labels] = load_mesh_gmsh_waterlevels(file_name)
% Load mesh data from a GMSH file with water levels.

% OUTPUTS:
%   coord: Coordinates of degrees of freedom, computed from mesh nodes and element_order.
%   elem: Element connectivity matrix (tetra2node), each row contains dof indices for an element.
%   surf: Surface connectivity matrix (triangle2node), each row contains dof indices for a surface.
%   Q: Logical arrat indicating dofs restricted by boundary conditions.
%   material: Material indices for each element

node = h5read(file_name, '/points');
tetra_cells = h5read(file_name, '/tetra_cells') + 1;  % Adjust for 1-based indexing in MATLAB
material = h5read(file_name, '/tetra_labels');
triangle_cells = h5read(file_name, '/triangles') + 1;  % Adjust for 1-based indexing in MATLAB
triangle_labels = h5read(file_name, '/triangle_labels');

Q = true(size(node));
% TODO: add stuff below, it is ok?

% Dirichlet bc u_x = 0 on x0, labels 7, 8

tmp = triangle_cells(:, triangle_labels == 7);
Q(1, tmp(:)) = false;  % x0 face
tmp = triangle_cells(:, triangle_labels == 8);
Q(1, tmp(:)) = false;  % x0 face

% Dirichlet bc u_x = 0 on x_max, label 9
tmp = triangle_cells(:, triangle_labels == 9);
Q(1, tmp(:)) = false;  % xmax fa

% Dirichlet bc u_y = 0 on y_0, label 10
tmp = triangle_cells(:, triangle_labels == 10);
Q(2, tmp(:)) = false;  % y0 face

% Dirichlet bc u_y = 0 on y_max, label 11
tmp = triangle_cells(:, triangle_labels == 11);
Q(2, tmp(:)) = false;  % ymax face

% Dirichlet bc u=0 on z0, label 12
tmp = triangle_cells(:, triangle_labels == 12);
Q(:, tmp(:)) = false;

% some random permutation of whatever is needed
coord = double(node([1 3 2], :));
Q = Q([1 3 2], :);

elem = double(tetra_cells);
elem=elem([1,2,3,4,5,6,7,10,9,8],:);
surf = double(triangle_cells);