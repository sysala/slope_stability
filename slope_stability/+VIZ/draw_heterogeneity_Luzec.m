function draw_heterogeneity_Luzec(coord,elem,mater)

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
  Q_black=(mater==0);
  Q_green=(mater==1);
  Q_cyan=(mater==2);
  Q_blue=(mater==3)|(mater==5)|(mater==7);
  Q_yellow=(mater==4);
  Q_magenta=(mater==6);
 
% visualization  
  figure
  hold on

  patch('Faces',elem(1:3,Q_black)','Vertices',coord',...
        'FaceColor','black','EdgeColor','none'); 
  patch('Faces',elem(1:3,Q_green)','Vertices',coord',...
        'FaceColor','green','EdgeColor','none'); 
  patch('Faces',elem(1:3,Q_cyan)','Vertices',coord',...
        'FaceColor','cyan','EdgeColor','none'); 
  patch('Faces',elem(1:3,Q_blue)','Vertices',coord',...
        'FaceColor','blue','EdgeColor','none');    
  patch('Faces',elem(1:3,Q_yellow)','Vertices',coord',...
        'FaceColor','yellow','EdgeColor','none'); 
  patch('Faces',elem(1:3,Q_magenta)','Vertices',coord',...
        'FaceColor','magenta','EdgeColor','none');    
  
  axis equal;  
%   axis off;
  hold off;

end