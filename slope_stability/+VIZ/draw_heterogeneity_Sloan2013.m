function draw_heterogeneity_Sloan2013(coord,elem,mater)

% =========================================================================
%
%  This function draws a given heterogeneity of the body
%
%  input data:
%    coord - coordinates of the nodes, size(coord)=(2,n_n) where n_n is a
%            number of nodes
%    elem - n_p x n_e array containing numbers of vertices defining each
%           element, n_e = number of elements
%    mater - 1 x n_e array containing numbers of materials in each
%           element, n_e = number of elements
%
% ======================================================================
%

% logical arrays specifying each material
  Q_blue=(mater==1);
  Q_cyan=(mater==0);
 
% visualization  
  figure
  hold on

  patch('Faces',elem(1:3,Q_blue)','Vertices',coord',...
        'FaceColor','blue','EdgeColor','none');    
  patch('Faces',elem(1:3,Q_cyan)','Vertices',coord',...
        'FaceColor','cyan','EdgeColor','none'); 
  
  axis equal;  
  axis off;
  hold off;

end