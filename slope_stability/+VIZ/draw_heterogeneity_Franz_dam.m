function draw_heterogeneity_Franz_dam(coord,elem,mater)

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
  Q_green=(mater==4); % Zone A = impervious core
  Q_yellow=(mater==3)|(mater==7); % Zone B = transition to dam fill
  Q_blue=(mater==8); % Zone C1 = earth fill downstream
  Q_cyan=(mater==2); % Zone C2 = earth fill upstream
  Q_red=(mater==9); % Zone D = drainage
  Q_magenta=(mater==6); % Zone E = grout curtain
  Q_black=(mater==1)|(mater==5); % Zone RF = rock fill 
  Q_grey=(mater==0); % Zone X = bedrock    
 
% visualization  
  figure
  hold on

  patch('Faces',elem(1:3,Q_green)','Vertices',coord',...
        'FaceColor','green','EdgeColor','none'); 
  patch('Faces',elem(1:3,Q_yellow)','Vertices',coord',...
        'FaceColor','yellow','EdgeColor','none'); 
  patch('Faces',elem(1:3,Q_blue)','Vertices',coord',...
        'FaceColor','blue','EdgeColor','none');    
  patch('Faces',elem(1:3,Q_cyan)','Vertices',coord',...
        'FaceColor','cyan','EdgeColor','none'); 
  patch('Faces',elem(1:3,Q_red)','Vertices',coord',...
        'FaceColor','red','EdgeColor','none'); 
  patch('Faces',elem(1:3,Q_magenta)','Vertices',coord',...
        'FaceColor','magenta','EdgeColor','none');
  patch('Faces',elem(1:3,Q_black)','Vertices',coord',...
        'FaceColor','black','EdgeColor','none');  
  patch('Faces',elem(1:3,Q_grey)','Vertices',coord',...
       'FaceColor',[1/2 1/2 1/2],'EdgeColor','none');    

  
  axis equal;  
%   axis off;
  hold off;

end