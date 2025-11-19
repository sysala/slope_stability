
function [coord,elem, mater]= mesh_P1_2D_Sloan2013...
                                (h,x1,x2,x3,y1,y2,y11,y12,y13,y21,y22,y23)

% =========================================================================
%
%  This function creates a triangular mesh for P1 elements and 
%  slope stability geometry
%
%  input data:
%    x1  - length of the body in front of the slope
%    x2  - length of the slope in x-direction
%    x3  - length of the body behind the slope
%    y1  - hight of the foundation
%    y2  - height of the slope
%    y11 - height of the foundation below the weaker zone
%    y12 - thickness of the weaker zone
%    y13 - height of the foundation above the weaker zone
%    y21 - heihgt of the water level next to the slope
%    y22 - difference between water levels on oposite slope sides
%    y23 - height of the slope above underground water level
%    h   - discretization parameter  
%
%  output data:
%    coord    - coordinates of the nodes, size(coord)=(2,n_n) where n_n is
%               a number of nodes
%    elem     - array containing numbers of nodes defining each element, 
%               size(elem)=(3,n_e), n_e = number of elements (triangles)
%
% ======================================================================
%

%
% Numbers of segments in x and y directions
%
  Nx12= round((x1+x2)/h); % slope segments + segments left from the slope
  Nx3 = round(x3/h) ;     % segments right from the slope 
  Nx = Nx12 + Nx3 ;       % total number of segments in x-direction
  Ny11 = ceil(y11/h) ;    
  Ny12 = ceil(y12/h) ;     
  Ny13 = ceil(y13/h) ;    
  Ny1 = Ny11+Ny12+Ny13 ;  % y-segments in the foundation
  Ny21 = ceil(y21/h) ;    
  Ny22 = ceil(y22/h) ;     
  Ny23 = ceil(y23/h) ;    
  Ny2 = Ny21+Ny22+Ny23 ;  % y-segments in the slope
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
  coord_y11=linspace(0,y11,Ny11+1);
  coord_y12=linspace(y11,y11+y12,Ny12+1);
  coord_y13=linspace(y11+y12,y1,Ny13+1);
  coord_y21=linspace(y1,y1+y21,Ny21+1);
  coord_y22=linspace(y1+y21,y1+y21+y22,Ny22+1);
  coord_y23=linspace(y1+y21+y22,y1+y2,Ny23+1);
  coord_y=[coord_y11 coord_y12(2:end) coord_y13(2:end) ...
           coord_y21(2:end) coord_y22(2:end) coord_y23(2:end)];
  
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
% The arrays elem and mater
%

  elem = zeros(3,2*Nx*Ny) ;
  mater = zeros(1,2*Nx*Ny) ;
  n_e = 0 ;
  
  % elements in foundation, below the weak layer
  for  j = 1:Ny11
    for  i = 1:Nx
      n_e = n_e + 1 ;
      elem(:,n_e) = [ V(i+1,j) ; V(i+1,j+1); V(i,j) ] ;
      n_e = n_e + 1 ;
      elem(:,n_e) = [ V(i,j+1) ; V(i,j)    ; V(i+1,j+1) ] ;
    end
  end

  % elements in the weak layer
  for  j = (Ny11+1):(Ny11+Ny12)
    for  i = 1:Nx
      n_e = n_e + 1 ;
      elem(:,n_e) = [ V(i+1,j) ; V(i+1,j+1); V(i,j) ] ;
      mater(n_e)=1;
      n_e = n_e + 1 ;
      elem(:,n_e) = [ V(i,j+1) ; V(i,j)    ; V(i+1,j+1) ] ;  
      mater(n_e)=1;      
    end
  end

  % elements in foundation, above the weak layer
  for  j = (Ny11+Ny12+1):Ny1
    for  i = 1:Nx
      n_e = n_e + 1 ;
      elem(:,n_e) = [ V(i+1,j) ; V(i+1,j+1); V(i,j) ] ;
      n_e = n_e + 1 ;
      elem(:,n_e) = [ V(i,j+1) ; V(i,j)    ; V(i+1,j+1) ] ;
    end
  end
   
  % elements left from the slope
  for  j = (Ny1+1):Ny
    for  i = 1:Nx12
      n_e = n_e + 1 ;
      elem(:,n_e) = [ V(i+1,j) ; V(i+1,j+1); V(i,j) ] ;
      n_e = n_e + 1 ;
      elem(:,n_e) = [ V(i,j+1) ; V(i,j)    ; V(i+1,j+1) ] ;
    end
  end

  %
  elem = elem(:,1:n_e) ;
  mater = mater(1:n_e) ;
  
end
