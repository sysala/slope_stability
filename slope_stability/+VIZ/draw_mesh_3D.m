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

  figure
  hold on
  patch('Faces',surf(1:3,:)','Vertices',coord','FaceVertexCData',...
        0*ones(size(coord,2),1),'FaceColor','white','EdgeColor','blue'); 
 
%   ind=unique(surf(:));
%   plot3( coord(1,ind),coord(2,ind),coord(3,ind), 'b.', 'MarkerSize',10);
  axis equal;  % real ratios
  view([0.5 1 -2]);
  hold off;
%  axis off;
VIZ.compute_bounding_edges(surf, coord, 1);


view([0.5 1 -2]);
box on;
axis equal;  % Maintain equal aspect ratios.
hold off;

% VIZ.compute_bounding_edges(surf, coord, 1);
% Set camera parameters for an orthographic projection.
ax2 = gca;
set(ax2,'Xdir','reverse');
cameraParams = struct( ...
    'Position', [1 4 -4] .* max(coord, [], 2)' * 2, ...
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