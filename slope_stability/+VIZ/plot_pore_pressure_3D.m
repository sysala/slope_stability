function fig = plot_pore_pressure_3D(pw, coord, surf)
%--------------------------------------------------------------------------
% plot_pore_pressure_3D plots the 3D pore pressure field on the mesh.
%
% This function uses the pore pressure field pw at each node,
% and then plots the field using the VIZ package.
%
% The function performs the following steps:
%   1. Shift the coordinates to ensure non-negative values.
%   2. Compute a scaling factor for visualization based on the pore pressure magnitude.
%   4. Call VIZ.draw_quantity_2D to plot the pore pressure field.
%
% INPUTS:
%   pw    - Pore pressure field (1 x n_n).
%   coord - Nodal coordinates (2 x n_n).
%   elem  - Element connectivity matrix.
%
% OUTPUT:
%   fig   - Handle to the generated figure.
%
%--------------------------------------------------------------------------
%%
% % Shift coordinates to ensure non-negative values.
% coord(1,:) = coord(1,:) - min(coord(1,:));
% coord(2,:) = coord(2,:) - min(coord(2,:));

% Plot displacements with color representing the norm of the displacement.
fig = VIZ.draw_quantity_3D(coord, surf, zeros(size(coord)), pw);
VIZ.compute_bounding_edges(surf, coord, 1);
title("Pore pressures [kPa]")
drawnow
end
