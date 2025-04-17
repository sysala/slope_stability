function fig = draw_quantity_2D(coord, surf, U, Q_node)
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


% -------------------------------------------------------------------------
% Plot all triangles at once using patch without edges.

if size(surf, 1) > 6
    % Local connectivity for the four sub-triangles in an element
    all_triangles = [1 6 4; 4 5 2; 3 5 6; 6 4 5];

    surf_all = [surf(all_triangles(1,:), :)';
        surf(all_triangles(2,:), :)';
        surf(all_triangles(3,:), :)';
        surf(all_triangles(4,:), :)'];


    patch('Faces', surf_all, 'Vertices', (coord + U)', ...
        'FaceVertexCData', Q_node(:), 'FaceColor', 'interp', 'EdgeColor', 'none');
elseif size(surf, 1) > 3
    % Local connectivity for the four sub-triangles in an element
    all_triangles = [1 6 5; 4 5 6; 2 4 6; 3 4 5];

    surf_all = [surf(all_triangles(1,:), :)';
        surf(all_triangles(2,:), :)';
        surf(all_triangles(3,:), :)';
        surf(all_triangles(4,:), :)'];


    patch('Faces', surf_all, 'Vertices', (coord + U)', ...
        'FaceVertexCData', Q_node(:), 'FaceColor', 'interp', 'EdgeColor', 'none');
else
    patch('Faces', surf', 'Vertices', (coord + U)', ...
        'FaceVertexCData', Q_node(:), 'FaceColor', 'interp', 'EdgeColor', 'none');

end
% Set axis properties
axis equal;
xlim([min(coord(1,:)), max(coord(1,:))]);
ylim([min(coord(2,:)), max(coord(2,:))]);
colorbar;
xlabel('X');    
ylabel('Y');
title("Deviatoric Strain");

drawnow;

end
