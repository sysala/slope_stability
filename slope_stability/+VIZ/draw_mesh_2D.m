function draw_mesh_2D(coord,elem)

% =========================================================================
%
%  This function draws a given finite element mesh 
%
%  input data:
%    coord - coordinates of the nodes, size(coord)=(2,n_n) where n_n is a
%            number of nodes
%    elem - n_p x n_e array containing numbers of nodes defining each
%           element, n_e = number of elements
%
% ======================================================================
%

  figure
  hold on

  patch('Faces',elem(1:3,:)','Vertices',coord','FaceColor','white',...
        'EdgeColor','black');     

%   plot( coord(1,:),coord(2,:), 'b.', 'MarkerSize',10);  
 
  axis equal;  
  axis off;
  hold off;

end