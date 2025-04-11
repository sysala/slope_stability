function [boundary_faces, coord_new, values_new] = boundary_stress_from_elements_3D(coord, elem, elem_values)
%--------------------------------------------------------------------------
% boundary_stress_from_elements_3D extracts the boundary faces and associated
% nodal values from a 3D finite element mesh.
%
% This function computes the boundary of a 3D mesh by assembling surface
% faces from the volumetric element connectivity and then determines which
% faces lie on the boundary. For each boundary face, the corresponding nodal 
% coordinates and element values (e.g., stress or strain values) are extracted.
%
% The procedure is as follows:
%   1. Four groups of faces are constructed from the element connectivity array
%      'elem' by selecting specific combinations of nodes.
%   2. Nodal coordinates (x, y, z) are extracted for each element.
%   3. The corresponding nodal values are similarly extracted from 'elem_values'.
%   4. A sparse mapping from faces to nodes is constructed, and boundary faces 
%      are identified as those which appear only once (i.e., their corresponding
%      face-to-face matrix has a maximum value of 3).
%   5. The boundary faces and their associated nodal coordinates and values are 
%      reshaped for further visualization.
%
% INPUTS:
%   coord       - Nodal coordinates, size (3, n_n), where n_n is the number of nodes.
%   elem        - Element connectivity matrix, where each column represents an element.
%   elem_values - Element-wise nodal values (e.g., stress or strain), arranged in
%                 the same order as the nodes in 'elem'.
%
% OUTPUTS:
%   boundary_faces - Matrix containing the connectivity of the boundary faces.
%   coord_new      - Coordinates of the nodes on the boundary, size (3, n_boundary_nodes).
%   values_new     - Nodal values associated with the boundary, as a row vector.
%
%--------------------------------------------------------------------------
%
% Construct four sets of faces from the element connectivity.
faces = [elem([1,2,3,5,6,7],:), elem([1,2,4,5,8,10],:), elem([1,3,4,7,9,10],:), elem([2,3,4,6,9,8],:)];

% Extract nodal coordinates for each element.
x_e = coord(1,:);
x_e = x_e(:);
x_e = x_e(elem);

y_e = coord(2,:);
y_e = y_e(:);
y_e = y_e(elem);

z_e = coord(3,:);
z_e = z_e(:);
z_e = z_e(elem);

% Extract corresponding element values.
tmp = elem_values;
values_f = [tmp([1,2,3,5,6,7],:), tmp([1,2,4,5,8,10],:), tmp([1,3,4,7,9,10],:), tmp([2,3,4,6,9,8],:)];

% Assemble nodal coordinate arrays for faces.
x_f = [x_e([1,2,3,5,6,7],:), x_e([1,2,4,5,8,10],:), x_e([1,3,4,7,9,10],:), x_e([2,3,4,6,9,8],:)];
y_f = [y_e([1,2,3,5,6,7],:), y_e([1,2,4,5,8,10],:), y_e([1,3,4,7,9,10],:), y_e([2,3,4,6,9,8],:)];
z_f = [z_e([1,2,3,5,6,7],:), z_e([1,2,4,5,8,10],:), z_e([1,3,4,7,9,10],:), z_e([2,3,4,6,9,8],:)];

% Create a sparse mapping from face indices to global node indices.
i_ind = repmat((1:size(faces,2)), 6, 1);
i_ind = i_ind(:);
j_ind = faces(:);
faces2coord = sparse(i_ind, j_ind, i_ind * 0 + 1);

% Construct a face-to-face matrix to identify boundary faces.
F2F = faces2coord * faces2coord';
F2F = F2F - diag(diag(F2F));

% Identify boundary faces: those with a maximum value of 3.
index_boundary_face = (max(F2F) == 3)';

% Extract boundary faces and associated nodal values.
boundary_faces = faces(:, index_boundary_face);
boundary_values = values_f(:, index_boundary_face);
boundary_x = x_f(:, index_boundary_face);
boundary_y = y_f(:, index_boundary_face);
boundary_z = z_f(:, index_boundary_face);

% For visualization purposes, reassign boundary_faces indices.
boundary_faces = reshape(1:numel(boundary_faces), 6, size(boundary_faces, 2));
coord_new = [boundary_x(:) boundary_y(:) boundary_z(:)]';
values_new = boundary_values(:)';

end
