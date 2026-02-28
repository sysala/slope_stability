function draw_mesh_3D(coord,surf)

% =========================================================================
%
%  This function draws mesh and nodal point on the surface of the body
%
%  input data:
%    coord - coordinates of the nodes, size(coord)=(3,n_n) where n_n is a
%            number of nodes
%    surf  - array containing numbers of nodes defining each surface element,
%            size(surf)=(n_p,n_s), n_s = number of surface elements
%
% ======================================================================
%

  fig = figure; %#ok<NASGU>
  hold on
  verts = coord';
  patch('Faces', surf(1:3,:)', 'Vertices', verts, ...
        'FaceColor', [0.93 0.93 0.93], 'FaceAlpha', 1.0, ...
        'EdgeColor', [0.30 0.45 0.80], 'LineWidth', 0.35);
 
%   ind=unique(surf(:));
%   plot3( coord(1,ind),coord(2,ind),coord(3,ind), 'b.', 'MarkerSize',10);
  axis equal;
  axis tight;
  ax2 = gca;
  VIZ.apply_standard_3d_camera(ax2, verts);
  hold off;
%  axis off;
VIZ.compute_bounding_edges(surf, coord, 1);


box off;
hold off;
drawnow;
end
