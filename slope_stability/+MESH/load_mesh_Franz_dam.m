function [coord, elem, Q, mater, surf] = load_mesh_Franz_dam(elem_type, path)
%LOAD_MESH_LUZEC Summary of this function goes here
%   Detailed explanation goes here
  
  % coordinates and Dirichlet nodes
  coord=load([path 'coordinates.txt'])';
  
  % arrays over finite elements
  elem=load([path 'elements.txt'])'+1; % vertices  
  mater=load([path 'materials.txt'])-1;  % material type
  % edges=load('edges.txt')'+1;            % vertices  
  % edge_types=load('edge_type.txt')'+1;   % vertices

  % detection of the surface nodes and extension of the arrays elem and 
  % coord in the case of P2 or P4 elements 
  if strcmp(elem_type,'P1')==1
    surf= MESH.surf_P1_2D(elem);
  end
  if strcmp(elem_type,'P2')==1
    [coord_mid, elem_mid, surf]= MESH.midpoints_P2(coord,elem);
    coord=[coord, coord_mid];
    elem=[elem; elem_mid];
  end
  if strcmp(elem_type,'P4')==1
    [coord_mid, elem_mid, surf]= MESH.midpoints_P4(coord,elem);
    coord=[coord, coord_mid];
    elem=[elem; elem_mid];
  end
  
  % logical array representing the Dirichlet boundary conditions
  x_max=max(coord(1,:));
  x_min=min(coord(1,:));
  y_min=min(coord(2,:));
  n_n=size(coord,2);
  Q=false(2,n_n);
  Q(1,:) = (coord(1,:)>x_min+0.2)&(coord(2,:)>y_min+0.2)&(coord(1,:)<x_max-0.2) ;
  Q(2,:) = coord(2,:)>y_min+0.2 ;

end

