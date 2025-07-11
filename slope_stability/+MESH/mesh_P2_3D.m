
function [coord,elem,surf,Q]=mesh_P2_3D(N_h,x1,x2,x3,y1,y2,z)

% =========================================================================
%
%  This function creates tetrahedral mesh for P2 elements
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
%              size(elem)=(10,n_e), n_e = number of elements
%    surf    - array containing numbers of nodes defining each surace element,
%              size(surf)=(6,n_s), n_s = number of surface elements
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
  n_node_xy = (2*N_x+1)*(2*N1_y+1)+(2*N2_x+1)*2*N2_y; % number of nodes in xy plane
  n_n = n_node_xy*(2*N_z+1);                % total number of nodes
  n_cell_xy = N_x*N1_y+N2_x*N2_y;           % number of cells in xy plane
  n_e = n_cell_xy*N_z*6;                    % total number of elements
   
%
% C - 3D auxilliary array that contains node numbers and that is important 
% for the mesh construction. Since the body is a union of two cuboids
% the array C also consists of two auxilliary 3D arrays, C1 and C2.
%
  C=zeros(2*N_x+1,2*N_y+1,2*N_z+1);
  C1=reshape(1:(2*N_x+1)*(2*N1_y+1)*(2*N_z+1),2*N_x+1,2*N1_y+1,2*N_z+1);
  C2=reshape(((2*N_x+1)*(2*N1_y+1)*(2*N_z+1)+1):n_n,2*N2_x+1,2*N2_y,2*N_z+1);
  C(1:(2*N_x+1),1:(2*N1_y+1),:)=C1;
  C((2*N1_x+1):(2*N_x+1),(2*N1_y+2):(2*N_y+1),:)=C2; 
  
%
% coordinates of nodes
%
  % coordinates in directions x, y and z
  coord_x=linspace(0,x1+x2+x3,2*N_x+1);
  coord_y=linspace(0,y1+y2,2*N_y+1);
  coord_z=linspace(0,z,2*N_z+1);
  %
  cy_x=x1+x2+x3+repmat(linspace(-x2-x3,0,2*N2_x+1)',1,2*N2_y).*...
       repmat(1-x2/((x2+x3)*y2)*coord_y(2:(2*N2_y+1)),2*N2_x+1,1);
 
  
  % long 1D arrays containing coordinates of all nodes in x,y,z directions
  c_x=x1+x2+x3-[repmat(coord_x,1,(2*N1_y+1)*(2*N_z+1)),...
         repmat(cy_x(:)',1,2*N_z+1) ];   
  c_y=[repmat(kron(coord_y(1:(2*N1_y+1)),ones(1,2*N_x+1)),1,2*N_z+1),...
           repmat(kron(coord_y((2*N1_y+2):(2*N_y+1)),ones(1,2*N2_x+1)),1,2*N_z+1)];     
  c_z=[kron(coord_z,ones(1,(2*N1_y+1)*(2*N_x+1))),...
           kron(coord_z,ones(1,(2*N2_x+1)*2*N2_y))];  
        
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

  V1=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V1(1:2:(2*N_x-1),1:2:(2*N1_y-1),1:2:(2*N_z-1))=1; 
  V1((2*N1_x+1):2:(2*N_x-1),(2*N1_y+1):2:(2*N_y-1) ,1:2:(2*N_z-1))=1;
  %
  V2=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V2(3:2:(2*N_x+1),1:2:(2*N1_y-1),1:2:(2*N_z-1))=1; 
  V2((2*N1_x+3):2:(2*N_x+1),(2*N1_y+1):2:(2*N_y-1),1:2:(2*N_z-1))=1;
  %
  V3=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V3(3:2:(2*N_x+1),3:2:(2*N1_y+1),1:2:(2*N_z-1))=1; 
  V3((2*N1_x+3):2:(2*N_x+1),(2*N1_y+3):2:(2*N_y+1),1:2:(2*N_z-1))=1;
  %
  V4=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V4(1:2:(2*N_x-1),3:2:(2*N1_y+1),1:2:(2*N_z-1))=1; 
  V4((2*N1_x+1):2:(2*N_x-1),(2*N1_y+3):2:(2*N_y+1),1:2:(2*N_z-1))=1;
  %
  V5=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V5(1:2:(2*N_x-1),1:2:(2*N1_y-1),3:2:(2*N_z+1))=1; 
  V5((2*N1_x+1):2:(2*N_x-1),(2*N1_y+1):2:(2*N_y-1),3:2:(2*N_z+1))=1;
  %
  V6=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V6(3:2:(2*N_x+1),1:2:(2*N1_y-1),3:2:(2*N_z+1))=1; 
  V6((2*N1_x+3):2:(2*N_x+1),(2*N1_y+1):2:(2*N_y-1),3:2:(2*N_z+1))=1;
  %
  V7=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V7(3:2:(2*N_x+1),3:2:(2*N1_y+1),3:2:(2*N_z+1))=1; 
  V7((2*N1_x+3):2:(2*N_x+1),(2*N1_y+3):2:(2*N_y+1),3:2:(2*N_z+1))=1;
  %
  V8=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V8(1:2:(2*N_x-1),3:2:(2*N1_y+1),3:2:(2*N_z+1))=1; 
  V8((2*N1_x+1):2:(2*N_x-1),(2*N1_y+3):2:(2*N_y+1),3:2:(2*N_z+1))=1;

  % logical arrays for midpoints, e.g. V12 represents the midpoints between
  % V1 and V2
  V12=false(2*N_x+1,2*N_y+1,2*N_z+1); 
  V12(2:2:(2*N_x),1:2:(2*N1_y-1),1:2:(2*N_z-1))=1; 
  V12((2*N1_x+2):2:(2*N_x),(2*N1_y+1):2:(2*N_y-1) ,1:2:(2*N_z-1))=1;  
  %
  V14=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V14(1:2:(2*N_x-1),2:2:(2*N1_y),1:2:(2*N_z-1))=1; 
  V14((2*N1_x+1):2:(2*N_x-1),(2*N1_y+2):2:(2*N_y) ,1:2:(2*N_z-1))=1;
  %
  V15=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V15(1:2:(2*N_x-1),1:2:(2*N1_y-1),2:2:(2*N_z))=1; 
  V15((2*N1_x+1):2:(2*N_x-1),(2*N1_y+1):2:(2*N_y-1),2:2:(2*N_z))=1;
  %
  V16=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V16(2:2:(2*N_x),1:2:(2*N1_y-1),2:2:(2*N_z))=1; 
  V16((2*N1_x+2):2:(2*N_x),(2*N1_y+1):2:(2*N_y-1) ,2:2:(2*N_z))=1;
  %
  V23=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V23(3:2:(2*N_x+1),2:2:(2*N1_y),1:2:(2*N_z-1))=1; 
  V23((2*N1_x+3):2:(2*N_x+1),(2*N1_y+2):2:(2*N_y),1:2:(2*N_z-1))=1;
  %
  V24=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V24(2:2:(2*N_x),2:2:(2*N1_y),1:2:(2*N_z-1))=1; 
  V24((2*N1_x+2):2:(2*N_x),(2*N1_y+2):2:(2*N_y),1:2:(2*N_z-1))=1;
  %
  V26=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V26(3:2:(2*N_x+1),1:2:(2*N1_y-1),2:2:(2*N_z))=1; 
  V26((2*N1_x+3):2:(2*N_x+1),(2*N1_y+1):2:(2*N_y-1),2:2:(2*N_z))=1;
  %
  V34=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V34(2:2:(2*N_x),3:2:(2*N1_y+1),1:2:(2*N_z-1))=1; 
  V34((2*N1_x+2):2:(2*N_x),(2*N1_y+3):2:(2*N_y+1),1:2:(2*N_z-1))=1;  
  %
  V36=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V36(3:2:(2*N_x+1),2:2:(2*N1_y),2:2:(2*N_z))=1; 
  V36((2*N1_x+3):2:(2*N_x+1),(2*N1_y+2):2:(2*N_y),2:2:(2*N_z))=1;
  %
  V37=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V37(3:2:(2*N_x+1),3:2:(2*N1_y+1),2:2:(2*N_z))=1; 
  V37((2*N1_x+3):2:(2*N_x+1),(2*N1_y+3):2:(2*N_y+1),2:2:(2*N_z))=1;
  %
  V45=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V45(1:2:(2*N_x-1),2:2:(2*N1_y),2:2:(2*N_z))=1; 
  V45((2*N1_x+1):2:(2*N_x-1),(2*N1_y+2):2:(2*N_y),2:2:(2*N_z))=1;
  %
  V46=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V46(2:2:(2*N_x),2:2:(2*N1_y),2:2:(2*N_z))=1; 
  V46((2*N1_x+2):2:(2*N_x),(2*N1_y+2):2:(2*N_y),2:2:(2*N_z))=1;
  %
  V47=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V47(2:2:(2*N_x),3:2:(2*N1_y+1),2:2:(2*N_z))=1; 
  V47((2*N1_x+2):2:(2*N_x),(2*N1_y+3):2:(2*N_y+1),2:2:(2*N_z))=1;
  %
  V48=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V48(1:2:(2*N_x-1),3:2:(2*N1_y+1),2:2:(2*N_z))=1; 
  V48((2*N1_x+1):2:(2*N_x-1),(2*N1_y+3):2:(2*N_y+1),2:2:(2*N_z))=1;
  %
  V56=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V56(2:2:(2*N_x),1:2:(2*N1_y-1),3:2:(2*N_z+1))=1; 
  V56((2*N1_x+2):2:(2*N_x),(2*N1_y+1):2:(2*N_y-1),3:2:(2*N_z+1))=1;
  %
  V58=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V58(1:2:(2*N_x-1),2:2:(2*N1_y),3:2:(2*N_z+1))=1; 
  V58((2*N1_x+1):2:(2*N_x-1),(2*N1_y+2):2:(2*N_y),3:2:(2*N_z+1))=1;
  %
  V67=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V67(3:2:(2*N_x+1),2:2:(2*N1_y),3:2:(2*N_z+1))=1; 
  V67((2*N1_x+3):2:(2*N_x+1),(2*N1_y+2):2:(2*N_y),3:2:(2*N_z+1))=1;
  %
  V68=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V68(2:2:(2*N_x),2:2:(2*N1_y),3:2:(2*N_z+1))=1; 
  V68((2*N1_x+2):2:(2*N_x),(2*N1_y+2):2:(2*N_y),3:2:(2*N_z+1))=1;
  %
  V78=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V78(2:2:(2*N_x),3:2:(2*N1_y+1),3:2:(2*N_z+1))=1; 
  V78((2*N1_x+2):2:(2*N_x),(2*N1_y+3):2:(2*N_y+1),3:2:(2*N_z+1))=1;

  % used division of the unit cube into 6 tetrahedrons:   
  %   V1 V2 V4 V6 V12 V24 V14 V26 V46 V16
  %   V1 V4 V5 V6 V14 V45 V15 V46 V56 V16
  %   V4 V5 V6 V8 V45 V56 V46 V58 V68 V48
  %   V2 V3 V4 V6 V23 V34 V24 V36 V46 V26
  %   V3 V6 V7 V4 V36 V67 V37 V46 V47 V34
  %   V4 V6 V7 V8 V46 V67 V47 V68 V78 V48
  % size(aux_elem)=(6*10,n_e/6)
  aux_elem=[C(V1)'; C(V2)'; C(V4)'; C(V6)'; C(V12)'; C(V24)'; C(V14)'; C(V26)'; C(V46)'; C(V16)';
            C(V1)'; C(V4)'; C(V5)'; C(V6)'; C(V14)'; C(V45)'; C(V15)'; C(V46)'; C(V56)'; C(V16)';
            C(V4)'; C(V5)'; C(V6)'; C(V8)'; C(V45)'; C(V56)'; C(V46)'; C(V58)'; C(V68)'; C(V48)';
            C(V2)'; C(V3)'; C(V4)'; C(V6)'; C(V23)'; C(V34)'; C(V24)'; C(V36)'; C(V46)'; C(V26)';
            C(V3)'; C(V6)'; C(V7)'; C(V4)'; C(V36)'; C(V67)'; C(V37)'; C(V46)'; C(V47)'; C(V34)';
            C(V4)'; C(V6)'; C(V7)'; C(V8)'; C(V46)'; C(V67)'; C(V47)'; C(V68)'; C(V78)'; C(V48)' ];
        
  % the array elem, size(elem)=(10,n_e)
  elem=reshape(aux_elem,10,n_e);     
  
%
% Surface of the body - the array "surf"
%
  
  % For each face of the body, we define the restriction C_s of the array C 
  % and logical 2D arrays V1_s,...,V24_s which enable to select appropriate
  % nodes from the array C_s. We consider the following ordering of the 
  % nodes within the unit square:
  %      V1_s -> [0 0], V2_s -> [1 0], V3_s -> [1 1], V4_s -> [0 1]
  % V12_s -> [1/2 0], V23_s -> [1 1/2], V34_s -> [1/2 1], V14_s -> [0 1/2], V13_s=V24_s -> [1/2 1/2]
  % Finally, we use the division of a rectangle into 2 triangles which is
  % in accordance to the division of a prism into 6 tetrahedrons, see above.

  % Face 1: y=0 (the bottom of the body)
  %   ordering of the nodes creating the unit square:
  %     V1_s -> [0 0 0], V2_s -> [1 0 0], V5_s -> [0 0 1], V6_s -> [1 0 1] 
  %   used division of the square into 2 triangles:   
  %     V2_s V1_s V6_s V16_s V26_s V12_s
  %     V5_s V6_s V1_s V16_s V15_s V56_s
  C_s=zeros(2*N_x+1,2*N_z+1);
  C_s(:,:)=C(:,1,:);  
  V1_s=false(2*N_x+1,2*N_z+1);  V1_s(1:2:(2*N_x-1),1:2:(2*N_z-1))=1; 
  V2_s=false(2*N_x+1,2*N_z+1);  V2_s(3:2:(2*N_x+1),1:2:(2*N_z-1))=1; 
  V5_s=false(2*N_x+1,2*N_z+1);  V5_s(1:2:(2*N_x-1),3:2:(2*N_z+1))=1; 
  V6_s=false(2*N_x+1,2*N_z+1);  V6_s(3:2:(2*N_x+1),3:2:(2*N_z+1))=1; 
  V12_s=false(2*N_x+1,2*N_z+1);  V12_s(2:2:(2*N_x),1:2:(2*N_z-1))=1;
  V15_s=false(2*N_x+1,2*N_z+1);  V15_s(1:2:(2*N_x-1),2:2:(2*N_z))=1;
  V16_s=false(2*N_x+1,2*N_z+1);  V16_s(2:2:(2*N_x),2:2:(2*N_z))=1;
  V26_s=false(2*N_x+1,2*N_z+1);  V26_s(3:2:(2*N_x+1),2:2:(2*N_z))=1; 
  V56_s=false(2*N_x+1,2*N_z+1);  V56_s(2:2:(2*N_x),3:2:(2*N_z+1))=1;
  aux_surf=[C_s(V2_s)'; C_s(V1_s)'; C_s(V6_s)'; C_s(V16_s)'; C_s(V26_s)'; C_s(V12_s)';
            C_s(V5_s)'; C_s(V6_s)'; C_s(V1_s)'; C_s(V16_s)'; C_s(V15_s)'; C_s(V56_s)' ];
  surf1=reshape(aux_surf,6,2*N_x*N_z);       

  % Face 2: y=y1 (the top of body in front of the slope)
  %   ordering of the nodes creating the unit square:
  %     V4_s -> [0 1 0], V3_s -> [1 1 0], V8_s -> [0 1 1], V7_s -> [1 1 1] 
  %   used division of the square into 2 triangles:   
  %     V3_s V7_s V4_s V47_s V34_s V37_s 
  %     V8_s V4_s V7_s V47_s V78_s V48_s
  C_s=zeros(2*N_x+1,2*N_z+1);
  C_s(:,:)=C(:,2*N1_y+1,:);  
  V4_s=false(2*N_x+1,2*N_z+1);  V4_s(1:2:(2*N1_x-1),1:2:(2*N_z-1))=1; 
  V3_s=false(2*N_x+1,2*N_z+1);  V3_s(3:2:(2*N1_x+1),1:2:(2*N_z-1))=1; 
  V8_s=false(2*N_x+1,2*N_z+1);  V8_s(1:2:(2*N1_x-1),3:2:(2*N_z+1))=1; 
  V7_s=false(2*N_x+1,2*N_z+1);  V7_s(3:2:(2*N1_x+1),3:2:(2*N_z+1))=1; 
  V34_s=false(2*N_x+1,2*N_z+1);  V34_s(2:2:(2*N1_x),1:2:(2*N_z-1))=1;
  V37_s=false(2*N_x+1,2*N_z+1);  V37_s(3:2:(2*N1_x+1),2:2:(2*N_z))=1;
  V47_s=false(2*N_x+1,2*N_z+1);  V47_s(2:2:(2*N1_x),2:2:(2*N_z))=1;
  V48_s=false(2*N_x+1,2*N_z+1);  V48_s(1:2:(2*N1_x-1),2:2:(2*N_z))=1;
  V78_s=false(2*N_x+1,2*N_z+1);  V78_s(2:2:(2*N1_x),3:2:(2*N_z+1))=1; 
  aux_surf=[C_s(V3_s)'; C_s(V7_s)'; C_s(V4_s)'; C_s(V47_s)'; C_s(V34_s)'; C_s(V37_s)';
            C_s(V8_s)'; C_s(V4_s)'; C_s(V7_s)'; C_s(V47_s)'; C_s(V78_s)'; C_s(V48_s)' ];
  surf2=reshape(aux_surf,6,2*N1_x*N_z);  

  % Face 3: y=y1+y2 (the top of the slope)
  %   ordering of the nodes creating the unit square:
  %     V4_s -> [0 1 0], V3_s -> [1 1 0], V8_s -> [0 1 1], V7_s -> [1 1 1] 
  %   used division of the square into 2 triangles:   
  %     V3_s V7_s V4_s V47_s V34_s V37_s 
  %     V8_s V4_s V7_s V47_s V78_s V48_s 
  C_s=zeros(2*N_x+1,2*N_z+1);
  C_s(:,:)=C(:,end,:);  
  V4_s=false(2*N_x+1,2*N_z+1);  V4_s((2*N1_x+1):2:(2*N_x-1),1:2:(2*N_z-1))=1; 
  V3_s=false(2*N_x+1,2*N_z+1);  V3_s((2*N1_x+3):2:(2*N_x+1),1:2:(2*N_z-1))=1; 
  V8_s=false(2*N_x+1,2*N_z+1);  V8_s((2*N1_x+1):2:(2*N_x-1),3:2:(2*N_z+1))=1; 
  V7_s=false(2*N_x+1,2*N_z+1);  V7_s((2*N1_x+3):2:(2*N_x+1),3:2:(2*N_z+1))=1; 
  V34_s=false(2*N_x+1,2*N_z+1);  V34_s((2*N1_x+2):2:(2*N_x),1:2:(2*N_z-1))=1;
  V37_s=false(2*N_x+1,2*N_z+1);  V37_s((2*N1_x+3):2:(2*N_x+1),2:2:(2*N_z))=1;
  V47_s=false(2*N_x+1,2*N_z+1);  V47_s((2*N1_x+2):2:(2*N_x),2:2:(2*N_z))=1;
  V48_s=false(2*N_x+1,2*N_z+1);  V48_s((2*N1_x+1):2:(2*N_x-1),2:2:(2*N_z))=1;
  V78_s=false(2*N_x+1,2*N_z+1);  V78_s((2*N1_x+2):2:(2*N_x),3:2:(2*N_z+1))=1; 
  aux_surf=[C_s(V3_s)'; C_s(V7_s)'; C_s(V4_s)'; C_s(V47_s)'; C_s(V34_s)'; C_s(V37_s)';
            C_s(V8_s)'; C_s(V4_s)'; C_s(V7_s)'; C_s(V47_s)'; C_s(V78_s)'; C_s(V48_s)' ];
  surf3=reshape(aux_surf,6,2*N2_x*N_z);  

  % Face 4: x=0 (the left face of body in front of the slope)
  %   ordering of the nodes creating the unit square:
  %     V1_s -> [0 0 0], V4_s -> [0 1 0], V5_s -> [0 0 1], V8_s -> [0 1 1] 
  %   used division of the square into 2 triangles:   
  %     V1_s V4_s V5_s V45_s V15_s V14_s 
  %     V8_s V5_s V4_s V45_s V48_s V58_s 
  C_s=zeros(2*N_y+1,2*N_z+1);
  C_s(:,:)=C(1,:,:);  
  V1_s=false(2*N_y+1,2*N_z+1);  V1_s(1:2:(2*N1_y-1),1:2:(2*N_z-1))=1; 
  V4_s=false(2*N_y+1,2*N_z+1);  V4_s(3:2:(2*N1_y+1),1:2:(2*N_z-1))=1; 
  V5_s=false(2*N_y+1,2*N_z+1);  V5_s(1:2:(2*N1_y-1),3:2:(2*N_z+1))=1; 
  V8_s=false(2*N_y+1,2*N_z+1);  V8_s(3:2:(2*N1_y+1),3:2:(2*N_z+1))=1; 
  V14_s=false(2*N_y+1,2*N_z+1);  V14_s(2:2:(2*N1_y),1:2:(2*N_z-1))=1; 
  V15_s=false(2*N_y+1,2*N_z+1);  V15_s(1:2:(2*N1_y-1),2:2:(2*N_z))=1; 
  V45_s=false(2*N_y+1,2*N_z+1);  V45_s(2:2:(2*N1_y),2:2:(2*N_z))=1; 
  V48_s=false(2*N_y+1,2*N_z+1);  V48_s(3:2:(2*N1_y+1),2:2:(2*N_z))=1; 
  V58_s=false(2*N_y+1,2*N_z+1);  V58_s(2:2:(2*N1_y),3:2:(2*N_z+1))=1; 
  aux_surf=[C_s(V1_s)'; C_s(V4_s)'; C_s(V5_s)'; C_s(V45_s)'; C_s(V15_s)'; C_s(V14_s)';
            C_s(V8_s)'; C_s(V5_s)'; C_s(V4_s)'; C_s(V45_s)'; C_s(V48_s)'; C_s(V58_s)' ];
  surf4=reshape(aux_surf,6,2*N1_y*N_z);     
   
  % Face 5: x=x1 (the face representing the slope)
  %   ordering of the nodes creating the unit square:
  %     V1_s -> [0 0 0], V4_s -> [0 1 0], V5_s -> [0 0 1], V8_s -> [0 1 1] 
  %   used division of the square into 2 triangles:   
  %     V1_s V4_s V5_s  V45_s V15_s V14_s 
  %     V8_s V5_s V4_s  V45_s V48_s V58_s 
  C_s=zeros(2*N_y+1,2*N_z+1);
  C_s(:,:)=C(2*N1_x+1,:,:);  
  V1_s=false(2*N_y+1,2*N_z+1);  V1_s((2*N1_y+1):2:(2*N_y-1),1:2:(2*N_z-1))=1; 
  V4_s=false(2*N_y+1,2*N_z+1);  V4_s((2*N1_y+3):2:(2*N_y+1),1:2:(2*N_z-1))=1; 
  V5_s=false(2*N_y+1,2*N_z+1);  V5_s((2*N1_y+1):2:(2*N_y-1),3:2:(2*N_z+1))=1; 
  V8_s=false(2*N_y+1,2*N_z+1);  V8_s((2*N1_y+3):2:(2*N_y+1),3:2:(2*N_z+1))=1; 
  V14_s=false(2*N_y+1,2*N_z+1);  V14_s((2*N1_y+2):2:(2*N_y),1:2:(2*N_z-1))=1; 
  V15_s=false(2*N_y+1,2*N_z+1);  V15_s((2*N1_y+1):2:(2*N_y-1),2:2:(2*N_z))=1; 
  V45_s=false(2*N_y+1,2*N_z+1);  V45_s((2*N1_y+2):2:(2*N_y),2:2:(2*N_z))=1; 
  V48_s=false(2*N_y+1,2*N_z+1);  V48_s((2*N1_y+3):2:(2*N_y+1),2:2:(2*N_z))=1; 
  V58_s=false(2*N_y+1,2*N_z+1);  V58_s((2*N1_y+2):2:(2*N_y),3:2:(2*N_z+1))=1; 
  aux_surf=[C_s(V1_s)'; C_s(V4_s)'; C_s(V5_s)'; C_s(V45_s)'; C_s(V15_s)'; C_s(V14_s)';
            C_s(V8_s)'; C_s(V5_s)'; C_s(V4_s)'; C_s(V45_s)'; C_s(V48_s)'; C_s(V58_s)' ];
  surf5=reshape(aux_surf,6,2*N2_y*N_z);  

  % Face 6: x=x1+x2+x3 (the face right from the slope)
  %   ordering of the nodes creating the unit square:
  %     V2_s -> [1 0 0], V3_s -> [1 1 0], V6_s -> [1 0 1], V7_s -> [1 1 1] 
  %   used division of the square into 2 triangles:   
  %     V2_s V6_s V3_s  V36_s V23_s V26_s 
  %     V7_s V3_s V6_s  V36_s V67_s V37_s 
  C_s=zeros(2*N_y+1,2*N_z+1);
  C_s(:,:)=C(end,:,:);  
  V2_s=false(2*N_y+1,2*N_z+1);  V2_s(1:2:(2*N_y-1),1:2:(2*N_z-1))=1; 
  V3_s=false(2*N_y+1,2*N_z+1);  V3_s(3:2:(2*N_y+1),1:2:(2*N_z-1))=1; 
  V6_s=false(2*N_y+1,2*N_z+1);  V6_s(1:2:(2*N_y-1),3:2:(2*N_z+1))=1; 
  V7_s=false(2*N_y+1,2*N_z+1);  V7_s(3:2:(2*N_y+1),3:2:(2*N_z+1))=1; 
  V23_s=false(2*N_y+1,2*N_z+1);  V23_s(2:2:(2*N_y),1:2:(2*N_z-1))=1;
  V26_s=false(2*N_y+1,2*N_z+1);  V26_s(1:2:(2*N_y-1),2:2:(2*N_z))=1;
  V36_s=false(2*N_y+1,2*N_z+1);  V36_s(2:2:(2*N_y),2:2:(2*N_z))=1;
  V37_s=false(2*N_y+1,2*N_z+1);  V37_s(3:2:(2*N_y+1),2:2:(2*N_z))=1;
  V67_s=false(2*N_y+1,2*N_z+1);  V67_s(2:2:(2*N_y),3:2:(2*N_z+1))=1;
  aux_surf=[C_s(V2_s)'; C_s(V6_s)'; C_s(V3_s)'; C_s(V36_s)'; C_s(V23_s)'; C_s(V26_s)';
            C_s(V7_s)'; C_s(V3_s)'; C_s(V6_s)'; C_s(V36_s)'; C_s(V67_s)'; C_s(V37_s)' ];
  surf6=reshape(aux_surf,6,2*N_y*N_z);  

  % Face 7: z=0 (the front face of the body)
  %   ordering of the nodes creating the unit square:
  %     V1_s -> [0 0 0], V2_s -> [1 0 0], V4_s -> [0 1 0], V3_s -> [1 1 0] 
  %   used division of the square into 2 triangles:   
  %     V1_s V2_s V4_s  V24_s V14_s V12_s 
  %     V3_s V4_s V2_s  V24_s V23_s V34_s 
  C_s=zeros(2*N_x+1,2*N_y+1);
  C_s(:,:)=C(:,:,1);  
  V1_s=false(2*N_x+1,2*N_y+1); V1_s(1:2:(2*N_x-1),1:2:(2*N1_y-1))=1; V1_s((2*N1_x+1):2:(2*N_x-1),(2*N1_y+1):2:(2*N_y-1))=1; 
  V2_s=false(2*N_x+1,2*N_y+1); V2_s(3:2:(2*N_x+1),1:2:(2*N1_y-1))=1; V2_s((2*N1_x+3):2:(2*N_x+1),(2*N1_y+1):2:(2*N_y-1))=1;
  V4_s=false(2*N_x+1,2*N_y+1); V4_s(1:2:(2*N_x-1),3:2:(2*N1_y+1))=1; V4_s((2*N1_x+1):2:(2*N_x-1),(2*N1_y+3):2:(2*N_y+1))=1;
  V3_s=false(2*N_x+1,2*N_y+1); V3_s(3:2:(2*N_x+1),3:2:(2*N1_y+1))=1; V3_s((2*N1_x+3):2:(2*N_x+1),(2*N1_y+3):2:(2*N_y+1))=1;
  V12_s=false(2*N_x+1,2*N_y+1); V12_s(2:2:(2*N_x),1:2:(2*N1_y-1))=1; V12_s((2*N1_x+2):2:(2*N_x),(2*N1_y+1):2:(2*N_y-1))=1; 
  V14_s=false(2*N_x+1,2*N_y+1); V14_s(1:2:(2*N_x-1),2:2:(2*N1_y))=1; V14_s((2*N1_x+1):2:(2*N_x-1),(2*N1_y+2):2:(2*N_y))=1; 
  V23_s=false(2*N_x+1,2*N_y+1); V23_s(3:2:(2*N_x+1),2:2:(2*N1_y))=1; V23_s((2*N1_x+3):2:(2*N_x+1),(2*N1_y+2):2:(2*N_y))=1;
  V24_s=false(2*N_x+1,2*N_y+1); V24_s(2:2:(2*N_x),2:2:(2*N1_y))=1; V24_s((2*N1_x+2):2:(2*N_x),(2*N1_y+2):2:(2*N_y))=1;
  V34_s=false(2*N_x+1,2*N_y+1); V34_s(2:2:(2*N_x),3:2:(2*N1_y+1))=1; V34_s((2*N1_x+2):2:(2*N_x),(2*N1_y+3):2:(2*N_y+1))=1;
  aux_surf=[C_s(V1_s)'; C_s(V2_s)'; C_s(V4_s)'; C_s(V24_s)'; C_s(V14_s)'; C_s(V12_s)';
            C_s(V3_s)'; C_s(V4_s)'; C_s(V2_s)'; C_s(V24_s)'; C_s(V23_s)'; C_s(V34_s)' ];
  surf7=reshape(aux_surf,6,2*N_x*N1_y+2*N2_x*N2_y);  
  
  % Face 8: z=z_max (the back face of the body)
  %   ordering of the nodes creating the unit square:
  %     V5_s -> [0 0 1], V6_s -> [1 0 1], V8_s -> [0 1 1], V7_s -> [1 1 1] 
  %   used division of the square into 2 triangles:   
  %     V5_s V8_s V6_s  V68_s V56_s V58_s 
  %     V7_s V6_s V8_s  V68_s V78_s V67_s 
  C_s=zeros(2*N_x+1,2*N_y+1);
  C_s(:,:)=C(:,:,end);  
  V5_s=false(2*N_x+1,2*N_y+1); V5_s(1:2:(2*N_x-1),1:2:(2*N1_y-1))=1; V5_s((2*N1_x+1):2:(2*N_x-1),(2*N1_y+1):2:(2*N_y-1))=1; 
  V6_s=false(2*N_x+1,2*N_y+1); V6_s(3:2:(2*N_x+1),1:2:(2*N1_y-1))=1; V6_s((2*N1_x+3):2:(2*N_x+1),(2*N1_y+1):2:(2*N_y-1))=1;
  V8_s=false(2*N_x+1,2*N_y+1); V8_s(1:2:(2*N_x-1),3:2:(2*N1_y+1))=1; V8_s((2*N1_x+1):2:(2*N_x-1),(2*N1_y+3):2:(2*N_y+1))=1;
  V7_s=false(2*N_x+1,2*N_y+1); V7_s(3:2:(2*N_x+1),3:2:(2*N1_y+1))=1; V7_s((2*N1_x+3):2:(2*N_x+1),(2*N1_y+3):2:(2*N_y+1))=1;
  V56_s=false(2*N_x+1,2*N_y+1); V56_s(2:2:(2*N_x),1:2:(2*N1_y-1))=1; V56_s((2*N1_x+2):2:(2*N_x),(2*N1_y+1):2:(2*N_y-1))=1; 
  V58_s=false(2*N_x+1,2*N_y+1); V58_s(1:2:(2*N_x-1),2:2:(2*N1_y))=1; V58_s((2*N1_x+1):2:(2*N_x-1),(2*N1_y+2):2:(2*N_y))=1; 
  V68_s=false(2*N_x+1,2*N_y+1); V68_s(2:2:(2*N_x),2:2:(2*N1_y))=1; V68_s((2*N1_x+2):2:(2*N_x),(2*N1_y+2):2:(2*N_y))=1;
  V67_s=false(2*N_x+1,2*N_y+1); V67_s(3:2:(2*N_x+1),2:2:(2*N1_y))=1; V67_s((2*N1_x+3):2:(2*N_x+1),(2*N1_y+2):2:(2*N_y))=1;
  V78_s=false(2*N_x+1,2*N_y+1); V78_s(2:2:(2*N_x),3:2:(2*N1_y+1))=1; V78_s((2*N1_x+2):2:(2*N_x),(2*N1_y+3):2:(2*N_y+1))=1;
  aux_surf=[C_s(V5_s)'; C_s(V8_s)'; C_s(V6_s)'; C_s(V68_s)'; C_s(V56_s)'; C_s(V58_s)';
            C_s(V7_s)'; C_s(V6_s)'; C_s(V8_s)'; C_s(V68_s)'; C_s(V78_s)'; C_s(V67_s)' ];
  surf8=reshape(aux_surf,6,2*N_x*N1_y+2*N2_x*N2_y);  
  
  % the array "surf"
  surf = [surf1 surf2 surf3 surf5 surf4 surf6 surf7 surf8] ;
  
  
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
