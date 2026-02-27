grho = 10;
file_path = 'meshes/comsol_mesh.h5';
[coord, elem, surf, Q, material_identifier, triangle_labels] = MESH.load_mesh_P2(file_path);
surf_nodes = unique(surf(:));
surf_points = coord(:, surf_nodes);
figure
plot3(surf_points(1,:), surf_points(2,:), surf_points(3,:), 'ro')
axis equal
[t,s] = title('Surf likely contains all boundaries, even those between domains with different parameters','Surf is a LIE!');
t.FontSize = 10;
s.FontSize = 16;

[Q_w, pw_D] = MESH.seepage_boundary_3D_hetero_comsol(coord, surf, triangle_labels, grho);
d_nodes = coord(:, Q_w == 0);
figure
plot3(d_nodes(1,:), d_nodes(2,:), d_nodes(3,:), 'ro')
title('Nodes on Dirichlet boundary')
axis equal

figure;
scatter3(coord(1,:), coord(2,:), coord(3,:), 36, pw_D, 'filled');
xlabel('X'); ylabel('Y'); zlabel('Z');
colorbar;
axis equal;
view(3);
colormap(jet);
title('Dirichlet boundary values')