function edges_merged = compute_bounding_edges(surf, coord, if_plot)

% surf: 3 x n matrix, each column contains indices of triangle nodes
% edges: m x 2 matrix, each row is a unique edge (sorted node indices)
surf_triangles = surf(1:3,:);


surf_triangles = sort(surf_triangles, 1);
surf_triangles = unique(surf_triangles', 'rows')';
% Get all edges from triangles
edges_raw = [surf_triangles([1 2],:), surf_triangles([2 3],:), surf_triangles([3 1],:)]; % 2 x (3*n)

% Sort node indices in each edge
edges_sorted = sort(edges_raw, 1);

% Transpose to get edges as rows
edges_list = edges_sorted';

% Get unique edges
edges = unique(edges_list, 'rows');

% Incidence matrix between triangles (surf) and nodes
num_tri = size(surf_triangles, 2);
num_nodes = max(surf_triangles(:));
tri_node_inc = sparse(reshape(repmat(1:num_tri, 3, 1), [], 1), surf_triangles(:), 1, num_tri, num_nodes);

% Incidence matrix between edges and nodes
num_edges = size(edges, 1);
edge_node_inc = sparse(reshape(repmat(1:num_edges, 1, 2), [], 1), edges(:), 1, num_edges, num_nodes);

% Incidence matrix between edges and triangles
edge_tri_inc = edge_node_inc * tri_node_inc';
edge_tri_inc = edge_tri_inc == 2; % Each edge is incident to a triangle if both nodes are in the triangle

[edge_tri,~,~] = find(edge_tri_inc');
edge_tri=reshape(edge_tri,2,[]);

% Compute unit normal for each triangle in surf_triangles
% Assume 'coord' is a 3 x n matrix of node coordinates

tri_coords1 = coord(:, surf_triangles(1,:))';
tri_coords2 = coord(:, surf_triangles(2,:))';
tri_coords3 = coord(:, surf_triangles(3,:))';

v1 = tri_coords2 - tri_coords1;
v2 = tri_coords3 - tri_coords1;

normals = cross(v1, v2, 2);
normals = normals ./ vecnorm(normals, 2, 2);

left_tri_norm = normals(edge_tri(1,:),:);
right_tri_norm = normals(edge_tri(2,:),:);

idx_edge = abs(sum(left_tri_norm.*right_tri_norm, 2))<(1-1e-3);

edges_on_edge = edges(idx_edge,:);

edges_merged = VIZ.merge_collinear_edges(edges_on_edge, coord, 3);

if if_plot
    hold on;
    num_edges_merged = size(edges_merged, 1);

    for i = 1:num_edges_merged
        pts = coord(:, edges_merged(i,:));
        plot3(pts(1,:), pts(2,:), pts(3,:), '-', 'Color', 'black', 'LineWidth', 1);
    end
end
end