function [coord_new, values_new] = refactor_stress_from_elements_2D(coord, elem, elem_values)
%--------------------------------------------------------------------------
% refactor_stress_from_elements_2D restructures element-wise stress or strain 
% values for visualization on a 2D finite element mesh.
%
% This function extracts nodal coordinates and associated element values 
% from a 2D triangular mesh. The element-wise values are rearranged for 
% visualization purposes, ensuring consistency in mapping stress or strain 
% data to individual elements.
%
% The procedure is as follows:
%   1. Extract nodal coordinates (x, y) for each element.
%   2. Reshape and organize these coordinates into a new format suitable 
%      for visualization.
%   3. Flatten the element-wise stress or strain values into a corresponding
%      row vector.
%
% INPUTS:
%   coord       - Nodal coordinates, size (2 x n_n), where n_n is the number of nodes.
%   elem        - Element connectivity matrix, where each column represents an element.
%   elem_values - Element-wise values (e.g., stress or strain), arranged in
%                 the same order as the nodes in 'elem'.
%
% OUTPUTS:
%   coord_new   - Reorganized nodal coordinates, size (2 x n_elem_nodes).
%   values_new  - Flattened element-wise values, stored as a row vector.
%
%--------------------------------------------------------------------------
%
% Extract nodal coordinates for each element.
x_f = coord(1, :);
x_f = x_f(:);
x_f = x_f(elem);

y_f = coord(2, :);
y_f = y_f(:);
y_f = y_f(elem);

% Rearrange coordinates into a 2D format.
coord_new = [x_f(:) y_f(:)]';

% Flatten element-wise values for visualization.
values_new = elem_values(:)';

end
