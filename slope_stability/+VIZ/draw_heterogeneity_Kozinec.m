function draw_heterogeneity_Kozinec(coord,elem,mater)

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
  Q_mater1=(mater==1);
  Q_mater2=(mater==2);
  Q_mater3=(mater==3);
  Q_mater4=(mater==4);
  Q_mater5=(mater==5);
  Q_mater6=(mater==6);
  Q_mater7=(mater==7);

% visualization  
  figure
  hold on

  patch('Faces',elem(1:3,Q_mater1)','Vertices',coord',...
        'FaceColor','cyan','EdgeColor','none'); 
  patch('Faces',elem(1:3,Q_mater2)','Vertices',coord',...
        'FaceColor','blue','EdgeColor','none'); 
  patch('Faces',elem(1:3,Q_mater3)','Vertices',coord',...
        'FaceColor','red','EdgeColor','none'); 
  patch('Faces',elem(1:3,Q_mater4)','Vertices',coord',...
        'FaceColor','black','EdgeColor','none');    
  patch('Faces',elem(1:3,Q_mater5)','Vertices',coord',...
        'FaceColor','red','EdgeColor','none'); 
  patch('Faces',elem(1:3,Q_mater6)','Vertices',coord',...
        'FaceColor','black','EdgeColor','none'); 
  patch('Faces',elem(1:3,Q_mater7)','Vertices',coord',...
        'FaceColor','yellow','EdgeColor','none');    
 
  axis equal;  
  axis off;
  hold off;

end