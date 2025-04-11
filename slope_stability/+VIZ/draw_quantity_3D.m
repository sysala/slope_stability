function fig = draw_quantity_3D(coord, surf, U, Q_node)
%--------------------------------------------------------------------------
% draw_quantity_3D visualizes a prescribed nodal quantity on a 3D mesh.
%
% This function depicts a nodal scalar quantity (Q_node) defined on a 3D mesh.
% It displays the deformed (or undeformed if U is zero) shape of the body by 
% adding the displacement U to the original coordinates. The mesh is provided 
% by the node coordinates (coord) and the surface connectivity (surf). The 
% nodal quantity Q_node is used to color the surface using interpolated colors.
%
% INPUTS:
%   coord   - Nodal coordinates, size (3, n_n), where n_n is the number of nodes.
%   surf    - Connectivity matrix for surface elements, size (n_p, n_s),
%             where n_p is the number of nodes per surface element and n_s is the 
%             number of surface elements.
%   U       - Nodal displacements, size (3, n_n). To visualize the undeformed shape,
%             set U = 0*U.
%   Q_node  - Prescribed nodal quantity to be visualized, size (1, n_n).
%
% OUTPUT:
%   fig     - Handle to the generated figure.
%
% NOTE:
%   The function adjusts the coordinates so that the minimum coordinate in each
%   direction is zero, then displays the mesh with interpolated coloring based on
%   Q_node. If the number of rows in 'surf' is greater than 3, the connectivity is
%   assumed to be for a higher-order element and is subdivided into triangles using
%   a fixed pattern.
%
%--------------------------------------------------------------------------
fig = figure;
hold on;

% Shift coordinates so that they are nonnegative.
coord(1,:) = coord(1,:) - min(coord(1,:));
coord(2,:) = coord(2,:) - min(coord(2,:));
coord(3,:) = coord(3,:) - min(coord(3,:));

% Visualization of the quantity.
if size(surf, 1) > 3
    % For higher order surface elements, subdivide into triangles.
    all_triangles = [1 4 6; 4 2 5; 4 5 6; 6 5 3];
    surf_all = [surf(all_triangles(1,:), :)'; 
                surf(all_triangles(2,:), :)'; 
                surf(all_triangles(3,:), :)'; 
                surf(all_triangles(4,:), :)'];
    patch('Faces', surf_all, 'Vertices', (coord + U)', ...
          'FaceVertexCData', Q_node', 'FaceColor', 'interp', 'EdgeColor', 'none');
else
    % For simple triangular surface elements.
    patch('Faces', surf(1:3, :)', 'Vertices', (coord + U)', ...
          'FaceVertexCData', Q_node', 'FaceColor', 'interp', 'EdgeColor', 'none');
end

colorbar;
view([0.5 1 -2]);
box on;
axis equal;  % Maintain equal aspect ratios.
hold off;

% Set camera parameters for an orthographic projection.
ax2 = gca;
cameraParams = struct( ...
    'Position', [-1 4 -4] .* max(coord, [], 2)' * 2, ...
    'Target', max(coord, [], 2)' / 2, ...
    'UpVector', [0 180 0], ...
    'ViewAngle', 0, ...
    'Projection', 'orthographic');

set(ax2, 'CameraPosition', cameraParams.Position);
set(ax2, 'CameraTarget', cameraParams.Target);
set(ax2, 'CameraUpVector', cameraParams.UpVector);
set(ax2, 'CameraViewAngle', cameraParams.ViewAngle);
set(ax2, 'Projection', cameraParams.Projection);

drawnow;

end
