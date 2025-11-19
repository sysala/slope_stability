function [coord, elem, surf, Q, material, triangle_labels] = load_mesh_gmsh_waterlevels(file_name)
% Load mesh data from a GMSH file with water levels.

% OUTPUTS:
%   coord: Coordinates of degrees of freedom, computed from mesh nodes and element_order.
%   elem: Element connectivity matrix (tetra2node), each row contains dof indices for an element.
%   surf: Surface connectivity matrix (triangle2node), each row contains dof indices for a surface.
%   Q: Logical array indicating dofs restricted by Dirichlet boundary conditions for displacements.
%   material: Material indices for each element.
%      0 - General foundation
%      1 - Relatively weak foundation
%      2 - General slope mass
%      3 - Cover layer
%   triangle_labels: An array specifying a surface type for each surface 
%                    triangle. The coordinate y is the height of the body.
%      1 - bottom layer
%      2 - middle layer
%      3 - top layer
%      4 - cover
%      5 - dry part of the slope face
%      6 - wet part of the slope face
%      7 - dry part of the face with x=0
%      8 - wet part of the face with x=0
%      9 - face with x=x_max
%     10 - face with z=0
%     11 - face with z=z_max
%     12 - face with y=0 - bottom of the foundation
%     13 - face with y=y_max - top of the slope
%     14 - water bed face (constant y representing height of the foundation) 

% reading the mesh data from the file
node = h5read(file_name, '/points');
tetra_cells = h5read(file_name, '/tetra_cells') + 1;  % Adjust for 1-based indexing in MATLAB
material = h5read(file_name, '/tetra_labels') - 1 ;   % Values 0, 1, 2, 3 
triangle_cells = h5read(file_name, '/triangles') + 1;  % Adjust for 1-based indexing in MATLAB
triangle_labels = h5read(file_name, '/triangle_labels');

% The array Q - initialization
Q = true(size(node));

% Dirichlet bc u_x = 0 on the face with x=0, labels 7, 8
tmp = triangle_cells(:, triangle_labels == 7);
Q(1, tmp(:)) = false;  
tmp = triangle_cells(:, triangle_labels == 8);
Q(1, tmp(:)) = false;  

% Dirichlet bc u_x = 0 on the face with x=x_max, label 9
tmp = triangle_cells(:, triangle_labels == 9);
Q(1, tmp(:)) = false;  

% Dirichlet bc u_z = 0 on the face with z=0, label 10
tmp = triangle_cells(:, triangle_labels == 10);
Q(2, tmp(:)) = false;  % 

% Dirichlet bc u_z = 0 on the face with z=z_max, label 11
tmp = triangle_cells(:, triangle_labels == 11);
Q(2, tmp(:)) = false;  

% Dirichlet bc u=0 at the bottom (y=0), label 12
tmp = triangle_cells(:, triangle_labels == 12);
Q(:, tmp(:)) = false;

% the arrays coord and Q
coord = double(node([1 3 2], :));
Q = Q([1 3 2], :);

% the array elem
elem = double(tetra_cells);
elem=elem([1,2,3,4,5,6,7,10,9,8],:);

% the array surf
surf = double(triangle_cells);