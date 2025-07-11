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
end