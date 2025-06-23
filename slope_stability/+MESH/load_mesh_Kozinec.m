function [coord, elem, Q, mater] = load_mesh_Kozinec(elem_type, path)
%LOAD_MESH_KOZINEC Summary of this function goes here
%   Detailed explanation goes here
  
  % coordinates and Dirichlet nodes
  coord=load([path 'coordinates3.txt'])';
  
  % arrays over finite elements
  elem=load([path 'elements3.txt'])'+1; % vertices  
  mater=load([path 'materials3.txt'])-1;  % material type

  % extensions of the arrays elem and coord in the case of P2 or P4 elements 
  if strcmp(elem_type,'P2')==1
    [coord_mid, elem_mid]= MESH.midpoints_P2(coord,elem);
    coord=[coord, coord_mid];
    elem=[elem; elem_mid];
  end
  if strcmp(elem_type,'P4')==1
    [coord_mid, elem_mid]= MESH.midpoints_P4(coord,elem);
    coord=[coord, coord_mid];
    elem=[elem; elem_mid];
  end
  
  % logical array representing the Dirichlet boundary conditions
  x_max=max(coord(1,:));
  x_min=min(coord(1,:));
  n_n=size(coord,2);
  Q=false(2,n_n);
  Q(1,:) = (coord(1,:)>x_min+0.2)&(coord(2,:)>0.2)&(coord(1,:)<x_max-0.2) ;
  Q(2,:) = coord(2,:)>0.2 ;
end

