
function [coord,elem,elem_ed,edge_el,Q]= mesh_P1_2D(h,x1,x2,x3,y1,y2)

% =========================================================================
%
%  This function creates a triangular mesh for P1 elements and 
%  slope stability geometry
%
%  input data:
%    x1 - length of the body in front of the slope
%    x2 - length of the slope in x-direction
%    x3 - length of the body behind the slope
%    y1 - hight of the body below the slope
%    y2 - height of the slope
%    h  - discretization parameter  
%
%  output data:
%    coord    - coordinates of the nodes, size(coord)=(2,n_n) where n_n is
%               a number of nodes
%    elem     - array containing numbers of nodes defining each element, 
%               size(elem)=(3,n_e), n_e = number of elements (triangles)
%    elem_ed  - array containing numbers of edges for each triangle. We
%               keep anticlockwise ordering of edges. In particular, if
%               V1, V2 and V3 are vertices of a triangle T then the
%               corresponding edges are denoted as E12, E23, E13. 
%               size(elem_edge)=(3,n_e)
%    edge_el  - array containing 2 numbers of elements adjacent to each
%               edge. We do not keep ordering of the triangles. The value 
%               zero means that the edge is on the boundary of the domain.
%               size(edge_elem)=(2,n_ed), n_ed = number of edges
%    Q        - logical array indicating the nodes where the homogeneous
%               Dirichlet boundary condition is considered, size(Q)=(2,n_n)
%
% ======================================================================
%

%
% Numbers of segments in x and y directions
%
  Nx12= round((x1+x2)/h); % slope segments + segments left from the slope
  Nx3 = round(x3/h) ;     % segments right from the slope 
  Nx = Nx12 + Nx3 ;       % total number of segments in x-direction
  Ny1 = round(y1/h) ;     % segments below the slope
  Ny2 = round(y2/h) ;     % slope segments 
  Ny = Ny1 + Ny2 ;        % total number of segments in y-direction

%
% The arrays coord and V. V is a 2D auxilliary array containing node 
% numbers used for  the construction of the array elem
%

  % coordinates in x-direction (below the slope)
  coord_x12=linspace(0,x1+x2,Nx12+1);
  coord_x3=linspace(x1+x2,x1+x2+x3,Nx3+1);
  coord_x=[coord_x12 coord_x3(2:end)];
  
  % coordinates in y-direction
  coord_y1=linspace(0,y1,Ny1+1);
  coord_y2=linspace(y1,y1+y2,Ny2+1);
  coord_y=[coord_y1 coord_y2(2:end)];
  
  % sizes of arrays coord and V
  coord=zeros(2,(Ny1+1)*(Nx+1)+Ny2*(Nx12+1));
  V=zeros(Nx+1,Ny+1); n_n=0;
  
  % assembling of coord and V below the slope
  for j=1:(Ny1+1)
    for i=1:(Nx+1)
         n_n = n_n+1 ;
         V(i,j) = n_n ;
         coord(:,n_n) = [coord_x(i); coord_y(j)] ;
    end
  end
  
  % assembling of coord and V left from the slope - we keep a constant
  % number of nodes in x-direction.
  for j=(Ny1+2):(Ny+1)
      x_max=x1+x2*(y1+y2-coord_y(j))/y2;
      coord_x=linspace(0,x_max,Nx12+1);
    for i=1:(Nx12+1)
         n_n = n_n+1 ;
         V(i,j) = n_n ;
         coord(:,n_n) = [coord_x(i); coord_y(j)] ;
    end
  end

%  
% The array elem
%

  elem = zeros(3,2*Nx*Ny) ;
  n_e = 0 ;
  
  % elements below the slope
  for  j = 1:Ny1
    for  i = 1:Nx
      n_e = n_e + 1 ;
      elem(:,n_e) = [ V(i,j)   ; V(i+1,j)   ; V(i,j+1) ] ;
      n_e = n_e + 1 ;
      elem(:,n_e) = [ V(i+1,j+1) ; V(i,j+1) ; V(i+1,j)   ] ;
    end
  end
  n1_e=n_e; % number of elements below the slope
  
  % elements left from the slope
  for  j = (Ny1+1):Ny
    for  i = 1:Nx12
      n_e = n_e + 1 ;
      elem(:,n_e) = [ V(i,j)   ; V(i+1,j)   ; V(i,j+1) ] ;
      n_e = n_e + 1 ;
      elem(:,n_e) = [ V(i+1,j+1) ; V(i,j+1) ; V(i+1,j)   ] ;
    end
  end
  elem = elem(:,1:n_e) ;
  n2_e=n_e-n1_e; % number of elements left from the slope
  
%
% The array elem_ed=[elem1_ed, elem2_ed], where elem1_ed and elem2_ed 
%           represent elements below and left from the slope, respectively. 
%

  % A) assembling of elem1_d:
  % Eh1, Ev1, Ed1 - 2D auxilliary arrays containing numbers of particular 
  % horizontal, vertical and diagonal edges used for the construction of
  % the array elem1_ed
  n_ed=0; n1_ed=Nx*(Ny1+1); 
  Eh1=reshape(1:n1_ed,Nx  ,Ny1+1)+n_ed;
  n_ed=n_ed+n1_ed; n1_ed=(Nx+1)*Ny1;
  Ev1=reshape(1:n1_ed,Nx+1,Ny1  )+n_ed;
  n_ed=n_ed+n1_ed; n1_ed=Nx*Ny1;
  Ed1=reshape(1:n1_ed,Nx  ,Ny1  )+n_ed;  
  %  E12,E23,E34,E14,E13 are 1 x (Nx*Ny) arrays containing numbers of edges 
  %  for each square cell of the triangulation. Ordering of the 
  %  edges within the unit square is the following:
  %  E12: [0 0]->[1 0], E23: [1 0]->[1 1], E34: [0 1]->[1 1]
  %  E14: [0 0]->[0 1], E24: [1 0]->[0 1]
  E12=Eh1(1: Nx   ,1: Ny1   ); E12=E12(:)';
  E23=Ev1(2:(Nx+1),1: Ny1   ); E23=E23(:)';
  E34=Eh1(1: Nx   ,2:(Ny1+1)); E34=E34(:)';
  E14=Ev1(1: Nx   ,1: Ny1   ); E14=E14(:)';
  E24=Ed1(:)'; 
  % used division of a square cell into 2 triangles:   
  %    E12;  E24;  E14;
  %    E34;  E24;  E23
  aux_elem_ed=[ E12;  E24;  E14;
                E34;  E24;  E23 ];
  elem1_ed=reshape(aux_elem_ed,3,n1_e);        
  
  % B) assembling of elem2_d:
  % Eh2, Ev2, Ed2 - 2D auxilliary arrays containing numbers of particular 
  % horizontal, vertical and diagonal edges used for the construction of
  % the array elem2_ed
  n_ed=n_ed+n1_ed; n1_ed=Nx12*Ny2; 
  Eh2=[Eh1(1:Nx12,end), reshape(1:n1_ed,Nx12,Ny2)+n_ed];
  n_ed=n_ed+n1_ed; n1_ed=(Nx12+1)*Ny2;
  Ev2=reshape(1:n1_ed,Nx12+1,Ny2)+n_ed;
  n_ed=n_ed+n1_ed; n1_ed=Nx12*Ny2;
  Ed2=reshape(1:n1_ed,Nx12  ,Ny2)+n_ed;  
  %  E12,E23,E34,E14,E13 are 1 x (Nx*Ny) arrays containing numbers of edges 
  %  for each square cell of the triangulation. Ordering of the 
  %  edges within the unit square is the following:
  %  E12: [0 0]->[1 0], E23: [1 0]->[1 1], E34: [0 1]->[1 1]
  %  E14: [0 0]->[0 1], E24: [1 0]->[0 1]
  E12=Eh2(1: Nx12   ,1: Ny2   ); E12=E12(:)';
  E23=Ev2(2:(Nx12+1),1: Ny2   ); E23=E23(:)';
  E34=Eh2(1: Nx12   ,2:(Ny2+1)); E34=E34(:)';
  E14=Ev2(1: Nx12   ,1: Ny2   ); E14=E14(:)';
  E24=Ed2(:)'; 
  % used division of a square cell into 2 triangles:   
  %    E12;  E24;  E14;
  %    E34;  E24;  E23
  aux_elem_ed=[ E12;  E24;  E14;
                E34;  E24;  E23 ];
  elem2_ed=reshape(aux_elem_ed,3,n2_e);      
  
  %
  elem_ed=[elem1_ed, elem2_ed];
  
%
% The array edge_el 
%

  % Tlb1,Trt1,Tlb2,Trt2 - 2D auxilliary arrays containing numbers of  
  % particular left-bottom and right-top triangles used for the 
  % construction of the array edge_el. Tlb1,Trt1 - triangles below the
  % slope, Tlb2,Trt2 - triangles left from the slope
  % The value 0 is used due to the presence of external edges.
  Tlb1=zeros(Nx+1,Ny1+1); Trt1=zeros(Nx+1,Ny1+1);
  Tlb1(1:Nx,1:Ny1)=reshape(1:2:(n1_e-1),Nx,Ny1);
  Tlb1(1:Nx12,end)=(1:2:(2*Nx12-1))+n1_e;
  Trt1(2:(Nx+1),2:(Ny1+1))=reshape(2:2:n1_e,Nx,Ny1);
  
  Tlb2=zeros(Nx12+1,Ny2+1); Trt2=zeros(Nx12+1,Ny2+1);
  Tlb2(1:Nx12,1:Ny2)=reshape((n1_e+1):2:(n_e-1),Nx12,Ny2);
  Trt2(2:(Nx12+1),2:(Ny2+1))=reshape((n1_e+2):2:n_e,Nx12,Ny2);
  Trt2(2:(Nx12+1),1)=(2:2:(2*Nx12))+n1_e-2*Nx;
   
  % edge1_el_h - a part of edge_el for horizontal edges below the slope. 
  %  T1h,T2h are 1 x Nx*(Ny1+1) arrays containing numbers of triangles for
  %  to each horizontal edge of the triangulation. Indices from T1h and T2h
  %  represent triangles above and below the edge, respectively.
  T1h=Tlb1(1:Nx,:); T1h=T1h(:)';
  T2h=Trt1(2:(Nx+1),:); T2h=T2h(:)';
  edge1_el_h=[T1h; T2h];
  
  % edge2_el_h - a part of edge_el for horizontal edges left from the slope. 
  %  T1h,T2h are 1 x Nx12*Ny2 arrays containing numbers of triangles for
  %  to each horizontal edge of the triangulation. Indices from T1h and T2h
  %  represent triangles above and below the edge, respectively.
  T1h=Tlb2(1:Nx12,2:(Ny2+1)); T1h=T1h(:)';
  T2h=Trt2(2:(Nx12+1),2:(Ny2+1)); T2h=T2h(:)';
  edge2_el_h=[T1h; T2h];
  
  % edge1_el_v - a part of edge_el for vertical edges below the slope.
  %  T1v,T2v are 1 x (Nx+1)*Ny1 arrays containing numbers of triangles for
  %  to each vertical edge of the triangulation. Indices from T1v and T2v
  %  represent left and right triangles from the edge, respectively.
  T1v=Trt1(:,2:(Ny1+1)); T1v=T1v(:)';
  T2v=Tlb1(:,1: Ny1   ); T2v=T2v(:)';
  edge1_el_v=[T1v; T2v];
  
  % edge2_el_v - a part of edge_el for vertical edges left from the slope.
  %  T1v,T2v are 1x (Nx12+1)*Ny2 arrays containing numbers of triangles for
  %  to each vertical edge of the triangulation. Indices from T1v and T2v
  %  represent left and right triangles from the edge, respectively.
  T1v=Trt2(:,2:(Ny2+1)); T1v=T1v(:)';
  T2v=Tlb2(:,1: Ny2   ); T2v=T2v(:)';
  edge2_el_v=[T1v; T2v];
  
  % edge1_el_d - a part of edge_el for diagonal edges below the slope.
  %  T1d,T2d are 1 x Nx*Ny1 arrays containing numbers of triangles for
  %  to each diagonal edge of the triangulation. Indices from T1d and T2d
  %  represent left and right triangles from the edge, respectively.
  T1d=Tlb1(1: Nx   ,1: Ny1   ); T1d=T1d(:)';
  T2d=Trt1(2:(Nx+1),2:(Ny1+1)); T2d=T2d(:)';
  edge1_el_d=[T1d; T2d];
  
  % edge2_el_d - a part of edge_el for diagonal edges left from the slope.
  %  T1d,T2d are 1 x Nx12*Ny2 arrays containing numbers of triangles for
  %  to each diagonal edge of the triangulation. Indices from T1d and T2d
  %  represent left and right triangles from the edge, respectively.
  T1d=Tlb2(1: Nx12   ,1: Ny2   ); T1d=T1d(:)';
  T2d=Trt2(2:(Nx12+1),2:(Ny2+1)); T2d=T2d(:)';
  edge2_el_d=[T1d; T2d];
  
  %
  edge_el=[edge1_el_h, edge1_el_v, edge1_el_d,...
           edge2_el_h, edge2_el_v, edge2_el_d] ;  

%
% Boundary conditions - the logical array Q
%   

  x_max=max(coord(1,:))-1e-9;
  Q=false(2,n_n);
  Q(1,:) = (coord(1,:)>0)&(coord(2,:)>0)&(coord(1,:)<x_max) ;
  Q(2,:) = coord(2,:)>0 ;
  
end
