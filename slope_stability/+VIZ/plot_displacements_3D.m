function fig = plot_displacements_3D(U, coord, elem, scale_factor)
%--------------------------------------------------------------------------
% plot_displacements_3D plots the 3D displacement field on the mesh.
%
% This function computes the norm of the displacement field U, associates 
% these values with elements, and then extracts the boundary faces and corresponding
% stress/displacement values for visualization using the VIZ package.
%
% The function performs the following steps:
%   1. Compute the displacement magnitude at each node.
%   2. Extract element-wise values from the nodal displacement magnitudes.
%   3. Compute the boundary faces, boundary node coordinates, and associated values
%      (used here to visualize displacement magnitude) via VIZ.boundary_stress_from_elements_3D.
%   4. Extract displacement components (x, y, z) on the element faces.
%   5. Normalize the coordinates by subtracting their minima.
%   6. Compute a scaling factor based on the maximum boundary value and the extent
%      of the mesh.
%   7. Call VIZ.draw_quantity_3D to plot the displacement field (scaled) on the mesh.
%
% INPUTS:
%   U            - Displacement field (3 x n_n), where each row corresponds to a spatial component.
%   coord        - Nodal coordinates (3 x n_n).
%   elem         - Element connectivity matrix.
%   scale_factor - (Optional) Scaling constant used in displacement scaling. Default is 0.05.
%
% OUTPUT:
%   fig   - Handle to the generated figure.
%
%--------------------------------------------------------------------------

if nargin < 4
    scale_factor = 0.05;
end

%%
U_size = sqrt(sum(U.^2));
elem_values = U_size(elem);
%%
[boundary_faces, coord_boundary, values_boundary] = VIZ.boundary_stress_from_elements_3D(coord, elem, elem_values);

U_tmp = U(1,:);
U_x_elem = U_tmp(elem);
[~, ~, U_x_face] = VIZ.boundary_stress_from_elements_3D(coord, elem, U_x_elem);

U_tmp = U(2,:);
U_y_elem = U_tmp(elem);
[~, ~, U_y_face] = VIZ.boundary_stress_from_elements_3D(coord, elem, U_y_elem);

U_tmp = U(3,:);
U_z_elem = U_tmp(elem);
[~, ~, U_z_face] = VIZ.boundary_stress_from_elements_3D(coord, elem, U_z_elem);
%%
% Shift coordinates to non-negative values.
coord(1,:) = coord(1,:) - min(coord(1,:));
coord(2,:) = coord(2,:) - min(coord(2,:));
coord(3,:) = coord(3,:) - min(coord(3,:));

% Compute a scaling factor based on the maximum boundary value and mesh size.
scale = max((max(coord(:)) * scale_factor) / values_boundary(:));

% Plot displacements with color representing the norm of the displacement.
fig = VIZ.draw_quantity_3D(coord_boundary, boundary_faces, [U_x_face; U_y_face; U_z_face] * scale, values_boundary);
title("Displacements (color = norm of displacement)")
end
