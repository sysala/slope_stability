
function eps=penalty_parameters_2D(coord,elem)
% =========================================================================
%
%  This function determines the penalization parameters for the unconfined 
%  seepage problem depending on the element size.
%
%  input data:
%    coord    - coordinates of the nodes, size(coord)=(2,n_n) where n_n is 
%               a number of nodes 
%    elem     - n_p x n_e array containing numbers of nodes defining each
%               element, n_e = number of elements. 
%
%  output data:
%    eps -  1 x n_e array containing the penalty parameters for each
%           element
%
% ======================================================================
% 

%
% Initialization:
%
  n_e=size(elem,2);   % number of elements (triangles)
  grho=9.81;          % specific weight of water
  eps=zeros(1,n_e);
  
%
% Procedure of refinement
%
  for i=1:n_e
      
      % vertices of the i-th element 
      W1 =elem(1,i); % \
      W2 =elem(2,i); % - vertices of Ti
      W3 =elem(3,i); % /
      % length of the element edges
      L12=sqrt((coord(1,W1)-coord(1,W2))^2+(coord(2,W1)-coord(2,W2))^2);
      L13=sqrt((coord(1,W1)-coord(1,W3))^2+(coord(2,W1)-coord(2,W3))^2);
      L23=sqrt((coord(1,W2)-coord(1,W3))^2+(coord(2,W2)-coord(2,W3))^2);
      eps(i)=grho*min([L12,L13,L23])/2;
      
  end % for cyclus
  
end
