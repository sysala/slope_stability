function fig = plot_deviatoric_strain_2D(U, coord, elem, B)
%--------------------------------------------------------------------------
% plot_deviatoric_strain_2D visualizes the deviatoric strain field on a 2D mesh.
%
% This function computes the element-wise deviatoric strain from the displacement
% field U and the strain-displacement matrix B using the function 
% VIZ.get_elem_stress_2D. The computed deviatoric strain is then visualized as 
% a color map over the mesh.
%
% INPUTS:
%   U     - Displacement field, size (2 x n_n), where n_n is the number of nodes.
%   coord - Nodal coordinates, size (2 x n_n).
%   elem  - Element connectivity matrix, size (3 x n), where n is the number of elements.
%   B     - Strain-displacement matrix used to compute the strain.
%
% OUTPUT:
%   fig   - Handle to the generated figure.
%
%--------------------------------------------------------------------------
%%
% Compute element-wise deviatoric strain values.
elem_values = VIZ.get_elem_stress_2D(U, B, elem);

% Reshape the element connectivity matrix for visualization.
if size(elem,1) == 15
    elem = elem(1:6,:);
    elem_new = [reshape(1:numel(elem), size(elem));zeros(9,size(elem,2))];
else
    elem_new = reshape(1:numel(elem), size(elem));
end
% Extract x-coordinates of element nodes and reshape for plotting.
x = coord(1, :);
x = x(:);
x_elem = x(elem);
x_elem = reshape(x_elem, 1, []);

% Extract y-coordinates of eleměent nodes and reshape for plotting.
y = coord(2, :);
y = y(:);
y_elem = y(elem);
y_elem = reshape(y_elem, 1, []);

% Create new coordinate matrix.
coord_new = [x_elem; y_elem];

% Flatten element strain values for visualization.
val_new = elem_values(:);

% Plot the deviatoric strain field.
fig = VIZ.draw_quantity_2D(coord_new, elem_new, coord_new * 0, val_new);
title("Deviatoric Strain");
drawnow
end
