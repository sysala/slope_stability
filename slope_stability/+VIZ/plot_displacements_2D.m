function fig = plot_displacements_2D(U, coord, elem)
%--------------------------------------------------------------------------
% plot_displacements_2D plots the 2D displacement field on the mesh.
%
% This function computes the norm of the displacement field U at each node,
% scales the displacements for visualization, and then plots the field 
% using the VIZ package.
%
% The function performs the following steps:
%   1. Compute the displacement magnitude at each node.
%   2. Shift the coordinates to ensure non-negative values.
%   3. Compute a scaling factor for visualization based on the displacement magnitude.
%   4. Call VIZ.draw_quantity_2D to plot the displacement field.
%
% INPUTS:
%   U     - Displacement field (2 x n_n), where each row corresponds to a spatial component.
%   coord - Nodal coordinates (2 x n_n).
%   elem  - Element connectivity matrix.
%
% OUTPUT:
%   fig   - Handle to the generated figure.
%
%--------------------------------------------------------------------------
%%
% Compute the displacement magnitude at each node.
U_size = sqrt(sum(U.^2));

%%
% Shift coordinates to ensure non-negative values.
coord(1,:) = coord(1,:) - min(coord(1,:));
coord(2,:) = coord(2,:) - min(coord(2,:));

% Compute a scaling factor for visualization.
scale = max(U_size(:)) / (max(coord(:)) * 0.05);

% Plot displacements with color representing the norm of the displacement.
fig = VIZ.draw_quantity_2D(coord, elem, U / scale, U_size);
title("Displacements (color = norm of displacement)")
drawnow
end
