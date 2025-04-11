function draw_mesh_2D(coord,elem,size_x,size_y)

% =========================================================================
%
%  This function draws a used mesh 
%
%  input data:
%    coord    - coordinates of the nodes, size(coord)=(2,n_n) where n_n is
%               a number of nodes
%    elem     - array containing numbers of nodes defining each element, 
%               size(elem)=(n_p,n_e), n_e = number of elements (triangles)
%
% ======================================================================
%

  figure
  hold on
  coord_aux = [coord;zeros(1,size(coord,2))];
  patch('Faces',elem(1:3,:)','Vertices',coord_aux','FaceVertexCData',...
           0*ones(size(coord,2),1),'FaceColor','white','EdgeColor','blue'); 
       
  %plot( coord(1,:),coord(2,:), 'b.', 'MarkerSize',10);
  axis([-0.1 size_x+0.1 -0.1 size_y+0.1])
  axis equal;  %realne pomery
  axis off;
  view(2);     %standartni pohled ve 2D
  hold off;

end