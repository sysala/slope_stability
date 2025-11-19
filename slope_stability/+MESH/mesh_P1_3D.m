
function [coord,elem,surf,Q]=mesh_P1_3D(N_h,x1,x2,x3,y1,y2,z)

% =========================================================================
%
%  This function creates tetrahedral mesh for P1 elements
%
%  input data (integers):
%    N_h      - an integer defining a density of a uniform mesh
%    x1       - length of the body in front of the slope
%    x2       - length of the the slope in x-direction
%    x3       - length of the body behind the slope
%    y1       - hight of the body below the slope
%    y2       - height of the slope
%    z        - length of the body in z-direction 
%    
% body before coordinate transformation:
%      (0,x1+x2+x3)x(0,y1)x(0,z)\cup (x1,x1+x2+x3)x(y1,y2)x(0,z)
%
%  output data:
%    coord   - coordinates of the nodes, size(coord)=(3,n_n) where n_n is a
%              number of nodes
%    elem    - array containing numbers of nodes defining each element,
%              size(elem)=(4,n_e), n_e = number of elements
%    surf    - array containing numbers of nodes defining each surace element,
%              size(surf)=(3,n_s), n_s = number of surface elements
%    Q       - logical array indicating the nodes where the homogeneous 
%              Dirichlet boundary condition is considered, size(Q)=(3,n_n)
%
% ======================================================================
%

%
% numbers of segments, nodes and elements
%

  N_x = (x1+x2+x3)*N_h; % number of segments in x direction
  N_y = (y1+y2)*N_h;    % number of segments in y direction
  N_z = z*N_h;          % number of segments in z direction
  N1_x=x1*N_h;          % number of segments in x direction in front of the slope    
  N2_x=N_x-N1_x;        % number of the remaining segments in x direction 
  N1_y=y1*N_h;          % number of segments in y direction below the slope
  N2_y=N_y-N1_y;        % number of segments in y direction below the slope
  % 
  n_node_xy = (N_x+1)*(N1_y+1)+(N2_x+1)*N2_y; % number of nodes in xy plane
  n_n = n_node_xy*(N_z+1);                % total number of nodes
  n_cell_xy = N_x*N1_y+N2_x*N2_y;         % number of cells in xy plane
  n_e = n_cell_xy*N_z*6;                  % total number of elements

%
% C - 3D auxilliary array that contains node numbers and that is important 
% for the mesh construction. Since the body is a union of two cuboids
% the array C also consists of two auxilliary 3D arrays, C1 and C2.
%
  C=zeros(N_x+1,N_y+1,N_z+1);
  C1=reshape(1:(N_x+1)*(N1_y+1)*(N_z+1),N_x+1,N1_y+1,N_z+1);
  C2=reshape(((N_x+1)*(N1_y+1)*(N_z+1)+1):n_n,N2_x+1,N2_y,N_z+1);
  C(1:(N_x+1),1:(N1_y+1),:)=C1;
  C((N1_x+1):(N_x+1),(N1_y+2):(N_y+1),:)=C2;
  
%
% coordinates of nodes
%
  % coordinates in directions x, y and z
  coord_x=linspace(0,x1+x2+x3,N_x+1);
  coord_y=linspace(0,y1+y2,N_y+1);
  coord_z=linspace(0,z,N_z+1);
  %
  cy_x=x1+x2+x3+repmat(linspace(-x2-x3,0,N2_x+1)',1,N2_y).*...
       repmat(1-x2/((x2+x3)*y2)*coord_y(2:N2_y+1),N2_x+1,1);
  
  % long 1D arrays containing coordinates of all nodes in x,y,z directions
  c_x=x1+x2+x3-[repmat(coord_x,1,(N1_y+1)*(N_z+1)),...
         repmat(cy_x(:)',1,N_z+1) ];   
  c_y=[repmat(kron(coord_y(1:(N1_y+1)),ones(1,N_x+1)),1,N_z+1),...
           repmat(kron(coord_y((N1_y+2):(N_y+1)),ones(1,N2_x+1)),1,N_z+1)];     
  c_z=[kron(coord_z,ones(1,(N1_y+1)*(N_x+1))),...
           kron(coord_z,ones(1,(N2_x+1)*N2_y))];  
 
         
  % the required array of coordinates, size(coord)=(3,n_n)   
  coord=[c_x; c_y; c_z] ;  

% 
% construction of the array elem
%
  % ordering of the nodes creating the unit cube:
  %  V1 -> [0 0 0], V2 -> [1 0 0], V3 -> [1 1 0], V4 -> [0 1 0]
  %  V5 -> [0 0 1], V6 -> [1 0 1], V7 -> [1 1 1], V8 -> [0 1 1]
  %  V1,...,V8 are logical 3D arrays which enable to select appropriate
  %  nodes from the array C.

  V1=false(N_x+1,N_y+1,N_z+1);
  V1(1:N_x,1:N1_y,1:N_z)=1; 
  V1((N1_x+1):N_x,(N1_y+1):N_y,1:N_z)=1;
  %
  V2=false(N_x+1,N_y+1,N_z+1);
  V2(2:(N_x+1),1:N1_y,1:N_z)=1; 
  V2((N1_x+2):(N_x+1),(N1_y+1):N_y,1:N_z)=1;
  %
  V3=false(N_x+1,N_y+1,N_z+1);
  V3(2:(N_x+1),2:(N1_y+1),1:N_z)=1; 
  V3((N1_x+2):(N_x+1),(N1_y+2):(N_y+1),1:N_z)=1;
  %
  V4=false(N_x+1,N_y+1,N_z+1);
  V4(1:N_x,2:(N1_y+1),1:N_z)=1; 
  V4((N1_x+1):N_x,(N1_y+2):(N_y+1),1:N_z)=1;
  %
  V5=false(N_x+1,N_y+1,N_z+1);
  V5(1:N_x,1:N1_y,2:(N_z+1))=1; 
  V5((N1_x+1):N_x,(N1_y+1):N_y,2:(N_z+1))=1;
  %
  V6=false(N_x+1,N_y+1,N_z+1);
  V6(2:(N_x+1),1:N1_y,2:(N_z+1))=1; 
  V6((N1_x+2):(N_x+1),(N1_y+1):N_y,2:(N_z+1))=1;
  %
  V7=false(N_x+1,N_y+1,N_z+1);
  V7(2:(N_x+1),2:(N1_y+1),2:(N_z+1))=1; 
  V7((N1_x+2):(N_x+1),(N1_y+2):(N_y+1),2:(N_z+1))=1;
  %
  V8=false(N_x+1,N_y+1,N_z+1);
  V8(1:N_x,2:(N1_y+1),2:(N_z+1))=1; 
  V8((N1_x+1):N_x,(N1_y+2):(N_y+1),2:(N_z+1))=1;
 
  % used division of a prism into 6 tetrahedrons:   
  %   V1 V2 V4 V6
  %   V1 V4 V5 V6
  %   V4 V5 V6 V8
  %   V2 V3 V4 V6
  %   V3 V6 V7 V4
  %   V4 V6 V7 V8
  % size(aux_elem)=(6*4,n_e/6)
  aux_elem=[C(V1)'; C(V2)'; C(V4)'; C(V6)';
            C(V1)'; C(V4)'; C(V5)'; C(V6)';
            C(V4)'; C(V5)'; C(V6)'; C(V8)';
            C(V2)'; C(V3)'; C(V4)'; C(V6)';
            C(V3)'; C(V6)'; C(V7)'; C(V4)';
            C(V4)'; C(V6)'; C(V7)'; C(V8)' ];
        
  % the array elem, size(elem)=(4,n_e)          
  elem=reshape(aux_elem,4,n_e);     

%
% Surface of the body - the array "surf"
%
  
  % For each face of the body, we define the restriction C_s of the array C 
  % and logical 2D arrays V1_s,...,V8_s which enable to select appropriate
  % nodes from the array C_s. We consider the following ordering of the 
  % nodes within the unit square:
  %   V1_s -> [0 0], V2_s -> [1 0], V3_s -> [1 1], V4_s -> [0 1]
  % Finally, we use the division of a rectangle into 2 triangles which is
  % in accordance to the division of a prism into 6 tetrahedrons, see above.

  % Face 1: y=0 (the bottom of the body)
  %   ordering of the nodes creating the unit square:
  %     V1_s -> [0 0 0], V2_s -> [1 0 0], V5_s -> [0 0 1], V6_s -> [1 0 1] 
  %   used division of the square into 2 triangles:   
  %     V2_s V1_s V6_s  and  V5_s V6_s V1_s 
  C_s=zeros(N_x+1,N_z+1);
  C_s(:,:)=C(:,1,:);  
  V1_s=false(N_x+1,N_z+1);  V1_s(1:N_x    ,1:N_z    )=1; 
  V2_s=false(N_x+1,N_z+1);  V2_s(2:(N_x+1),1:N_z    )=1; 
  V5_s=false(N_x+1,N_z+1);  V5_s(1:N_x    ,2:(N_z+1))=1; 
  V6_s=false(N_x+1,N_z+1);  V6_s(2:(N_x+1),2:(N_z+1))=1; 
  aux_surf=[C_s(V2_s)'; C_s(V1_s)'; C_s(V6_s)';
            C_s(V5_s)'; C_s(V6_s)'; C_s(V1_s)' ];
  surf1=reshape(aux_surf,3,2*N_x*N_z);       

  % Face 2: y=y1 (the top of body in front of the slope)
  %   ordering of the nodes creating the unit square:
  %     V4_s -> [0 1 0], V3_s -> [1 1 0], V8_s -> [0 1 1], V7_s -> [1 1 1] 
  %   used division of the square into 2 triangles:   
  %     V3_s V7_s V4_s  and  V8_s V4_s V7_s 
  C_s=zeros(N_x+1,N_z+1);
  C_s(:,:)=C(:,N1_y+1,:);  
  V4_s=false(N_x+1,N_z+1);  V4_s(1:N1_x    ,1:N_z    )=1; 
  V3_s=false(N_x+1,N_z+1);  V3_s(2:(N1_x+1),1:N_z    )=1; 
  V8_s=false(N_x+1,N_z+1);  V8_s(1:N1_x    ,2:(N_z+1))=1; 
  V7_s=false(N_x+1,N_z+1);  V7_s(2:(N1_x+1),2:(N_z+1))=1; 
  aux_surf=[C_s(V3_s)'; C_s(V7_s)'; C_s(V4_s)';
            C_s(V8_s)'; C_s(V4_s)'; C_s(V7_s)' ];
  surf2=reshape(aux_surf,3,2*N1_x*N_z);  

  % Face 3: y=y1+y2 (the top of the slope)
  %   ordering of the nodes creating the unit square:
  %     V4_s -> [0 1 0], V3_s -> [1 1 0], V8_s -> [0 1 1], V7_s -> [1 1 1] 
  %   used division of the square into 2 triangles:   
  %     V3_s V7_s V4_s  and  V8_s V4_s V7_s 
  C_s=zeros(N_x+1,N_z+1);
  C_s(:,:)=C(:,end,:);  
  V4_s=false(N_x+1,N_z+1);  V4_s((N1_x+1):N_x    ,1:N_z    )=1; 
  V3_s=false(N_x+1,N_z+1);  V3_s((N1_x+2):(N_x+1),1:N_z    )=1; 
  V8_s=false(N_x+1,N_z+1);  V8_s((N1_x+1):N_x    ,2:(N_z+1))=1; 
  V7_s=false(N_x+1,N_z+1);  V7_s((N1_x+2):(N_x+1),2:(N_z+1))=1; 
  aux_surf=[C_s(V3_s)'; C_s(V7_s)'; C_s(V4_s)';
            C_s(V8_s)'; C_s(V4_s)'; C_s(V7_s)' ];
  surf3=reshape(aux_surf,3,2*N2_x*N_z);  

  % Face 4: x=0 (the left face of body in front of the slope)
  %   ordering of the nodes creating the unit square:
  %     V1_s -> [0 0 0], V4_s -> [0 1 0], V5_s -> [0 0 1], V8_s -> [0 1 1] 
  %   used division of the square into 2 triangles:   
  %     V1_s V4_s V5_s  and  V8_s V5_s V4_s 
  C_s=zeros(N_y+1,N_z+1);
  C_s(:,:)=C(1,:,:);  
  V1_s=false(N_y+1,N_z+1);  V1_s(1:N1_y    ,1:N_z    )=1; 
  V4_s=false(N_y+1,N_z+1);  V4_s(2:(N1_y+1),1:N_z    )=1; 
  V5_s=false(N_y+1,N_z+1);  V5_s(1:N1_y    ,2:(N_z+1))=1; 
  V8_s=false(N_y+1,N_z+1);  V8_s(2:(N1_y+1),2:(N_z+1))=1; 
  aux_surf=[C_s(V1_s)'; C_s(V4_s)'; C_s(V5_s)';
            C_s(V8_s)'; C_s(V5_s)'; C_s(V4_s)' ];
  surf4=reshape(aux_surf,3,2*N1_y*N_z);     
   
  % Face 5: x=x1 (the face representing the slope)
  %   ordering of the nodes creating the unit square:
  %     V1_s -> [0 0 0], V4_s -> [0 1 0], V5_s -> [0 0 1], V8_s -> [0 1 1] 
  %   used division of the square into 2 triangles:   
  %     V1_s V4_s V5_s  and  V8_s V5_s V4_s 
  C_s=zeros(N_y+1,N_z+1);
  C_s(:,:)=C(N1_x+1,:,:);  
  V1_s=false(N_y+1,N_z+1);  V1_s((N1_y+1):N_y    ,1:N_z    )=1; 
  V4_s=false(N_y+1,N_z+1);  V4_s((N1_y+2):(N_y+1),1:N_z    )=1; 
  V5_s=false(N_y+1,N_z+1);  V5_s((N1_y+1):N_y    ,2:(N_z+1))=1; 
  V8_s=false(N_y+1,N_z+1);  V8_s((N1_y+2):(N_y+1),2:(N_z+1))=1; 
  aux_surf=[C_s(V1_s)'; C_s(V4_s)'; C_s(V5_s)';
            C_s(V8_s)'; C_s(V5_s)'; C_s(V4_s)' ];
  surf5=reshape(aux_surf,3,2*N2_y*N_z);  

  % Face 6: x=x1+x2+x3 (the face right from the slope)
  %   ordering of the nodes creating the unit square:
  %     V2_s -> [1 0 0], V3_s -> [1 1 0], V6_s -> [1 0 1], V7_s -> [1 1 1] 
  %   used division of the square into 2 triangles:   
  %     V2_s V6_s V3_s  and  V7_s V3_s V6_s 
  C_s=zeros(N_y+1,N_z+1);
  C_s(:,:)=C(end,:,:);  
  V2_s=false(N_y+1,N_z+1);  V2_s(1:N_y    ,1:N_z    )=1; 
  V3_s=false(N_y+1,N_z+1);  V3_s(2:(N_y+1),1:N_z    )=1; 
  V6_s=false(N_y+1,N_z+1);  V6_s(1:N_y    ,2:(N_z+1))=1; 
  V7_s=false(N_y+1,N_z+1);  V7_s(2:(N_y+1),2:(N_z+1))=1; 
  aux_surf=[C_s(V2_s)'; C_s(V6_s)'; C_s(V3_s)';
            C_s(V7_s)'; C_s(V3_s)'; C_s(V6_s)' ];
  surf6=reshape(aux_surf,3,2*N_y*N_z);  

  % Face 7: z=0 (the front face of the body)
  %   ordering of the nodes creating the unit square:
  %     V1_s -> [0 0 0], V2_s -> [1 0 0], V4_s -> [0 1 0], V3_s -> [1 1 0] 
  %   used division of the square into 2 triangles:   
  %     V1_s V2_s V4_s  and  V3_s V4_s V2_s 
  C_s=zeros(N_x+1,N_y+1);
  C_s(:,:)=C(:,:,1);  
  V1_s=false(N_x+1,N_y+1); V1_s(1:N_x    ,1:N1_y    )=1; V1_s((N1_x+1):N_x    ,(N1_y+1):N_y    )=1; 
  V2_s=false(N_x+1,N_y+1); V2_s(2:(N_x+1),1:N1_y    )=1; V2_s((N1_x+2):(N_x+1),(N1_y+1):N_y)=1;
  V4_s=false(N_x+1,N_y+1); V4_s(1:N_x    ,2:(N1_y+1))=1; V4_s((N1_x+1):N_x    ,(N1_y+2):(N_y+1))=1;
  V3_s=false(N_x+1,N_y+1); V3_s(2:(N_x+1),2:(N1_y+1))=1; V3_s((N1_x+2):(N_x+1),(N1_y+2):(N_y+1))=1;
  aux_surf=[C_s(V1_s)'; C_s(V2_s)'; C_s(V4_s)';
            C_s(V3_s)'; C_s(V4_s)'; C_s(V2_s)' ];
  surf7=reshape(aux_surf,3,2*N_x*N1_y+2*N2_x*N2_y);  
  
  % Face 8: z=z_max (the back face of the body)
  %   ordering of the nodes creating the unit square:
  %     V5_s -> [0 0 1], V6_s -> [1 0 1], V8_s -> [0 1 1], V7_s -> [1 1 1] 
  %   used division of the square into 2 triangles:   
  %     V5_s V8_s V6_s  and  V7_s V6_s V8_s 
  C_s=zeros(N_x+1,N_y+1);
  C_s(:,:)=C(:,:,end);  
  V5_s=false(N_x+1,N_y+1); V5_s(1:N_x    ,1:N1_y    )=1; V5_s((N1_x+1):N_x    ,(N1_y+1):N_y    )=1; 
  V6_s=false(N_x+1,N_y+1); V6_s(2:(N_x+1),1:N1_y    )=1; V6_s((N1_x+2):(N_x+1),(N1_y+1):N_y)=1;
  V8_s=false(N_x+1,N_y+1); V8_s(1:N_x    ,2:(N1_y+1))=1; V8_s((N1_x+1):N_x    ,(N1_y+2):(N_y+1))=1;
  V7_s=false(N_x+1,N_y+1); V7_s(2:(N_x+1),2:(N1_y+1))=1; V7_s((N1_x+2):(N_x+1),(N1_y+2):(N_y+1))=1;
  aux_surf=[C_s(V5_s)'; C_s(V8_s)'; C_s(V6_s)';
            C_s(V7_s)'; C_s(V6_s)'; C_s(V8_s)' ];
  surf8=reshape(aux_surf,3,2*N_x*N1_y+2*N2_x*N2_y);  
  
  % the array "surf"
  surf = [surf1 surf2 surf3 surf4 surf5 surf6 surf7 surf8] ;
 
%
% Dirichlet boundary conditions
%
  
  % logical array indicating the nodes with the Dirichlet boundary cond.
%   Q = coord>0 ;
%   Q(1,(coord(1,:)==x1+x2+x3)) = 0;     
%   Q(3,(coord(3,:)==z)) = 0;     

  Q=false(3,n_n);
  Q(1,:) = (coord(1,:)>0)&(coord(2,:)>0)&(coord(1,:)<x1+x2+x3) ;
  Q(2,:) = coord(2,:)>0 ;
  Q(3,:) = (coord(3,:)>0)&(coord(2,:)>0)&(coord(3,:)<z) ;
end
