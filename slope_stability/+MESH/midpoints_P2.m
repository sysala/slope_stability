function [coord_mid, elem_mid]= midpoints_P2(coord,elem)

%  The aim of this function is to find numbers of midpoints and their
%  coordinates. It enables us to use the mesh for P2 elements.
%
%  input data:
%    coord - coordinates of the vertices, size(coord)=(2,n_n) where n_n is 
%            a number of vertices
%    elem - 3 x n_e array containing numbers of nodes defining each
%           element, n_e = number of elements. We keep anticlockwise
%           ordering of the nodes creating a particular element. 
%
%  output data:
%    coord_mid - coordinates of midpoints, size(coord)=(2,n_edge)
%    elem_mid - 3 x n_e array containing numbers of midpoints in each
%               element
%
% ======================================================================
%
  
  % numbers of elements and vertices
  n_e=size(elem,2);
  n_n=size(coord,2);
      
  % predefinition of unknown arrays
  coord_mid=zeros(2,2*n_e);
  elem_mid=zeros(3,n_e);
  
  % for cyclus over elements
  ind=0; % enlarging index specifying midpoints
  for i=1:n_e
      
      % vertices defining the i-th element 
      V1=elem(1,i);
      V2=elem(2,i);
      V3=elem(3,i);
      
      % analysis of the edge V2-V3
      if elem_mid(1,i)==0
       % creation of a new midpoint   
       ind=ind+1; % number of a new midpoint
       coord_mid(:,ind)=(coord(:,V2)+coord(:,V3))/2; % its coordinates
       elem_mid(1,i)=n_n+ind;   
       % finding the adjacent element j to i which contains the edge V2-V3
       [row1,col1]=find(elem==V2);
       [row2,col2]=find(elem==V3);
       j=setdiff(intersect(col1,col2),i);
       if ~isempty(j)
          % This case means that the edge V2-V3 is the intersection of the
          % elements i and j.
          v=elem(:,j);
          if V3==v(1)
              elem_mid(3,j)=n_n+ind;
          elseif V3==v(2)
              elem_mid(1,j)=n_n+ind;
          else
              elem_mid(2,j)=n_n+ind;
          end
       end
      end % if (for the edge V2-V3)
      
      % analysis of the edge V3-V1
      if elem_mid(2,i)==0
       % creation of a new midpoint   
       ind=ind+1; % number of a new midpoint
       coord_mid(:,ind)=(coord(:,V3)+coord(:,V1))/2; % its coordinates
       elem_mid(2,i)=n_n+ind;             
       % finding the adjacent element j to i which contains the edge V3-V1
       [row1,col1]=find(elem==V3);
       [row2,col2]=find(elem==V1);
       j=setdiff(intersect(col1,col2),i);
       if ~isempty(j)
          % This case means that the edge V2-V3 is the intersection of the
          % elements i and j.
          v=elem(:,j);
          if V1==v(1)
              elem_mid(3,j)=n_n+ind;
          elseif V1==v(2)
              elem_mid(1,j)=n_n+ind;
          else
              elem_mid(2,j)=n_n+ind;
          end
       end
      end % if (for the edge V3-V1)
      
      % analysis of the edge V1-V2
      if elem_mid(3,i)==0
       % creation of a new midpoint   
       ind=ind+1; % number of a new midpoint
       coord_mid(:,ind)=(coord(:,V1)+coord(:,V2))/2; % its coordinates
       elem_mid(3,i)=n_n+ind;   
       % finding the adjacent element j to i which contains the edge V1-V2
       [row1,col1]=find(elem==V1);
       [row2,col2]=find(elem==V2);
       j=setdiff(intersect(col1,col2),i);
       if ~isempty(j)
          % This case means that the edge V2-V3 is the intersection of the
          % elements i and j.
          v=elem(:,j);
          if V2==v(1)
              elem_mid(3,j)=n_n+ind;
          elseif V2==v(2)
              elem_mid(1,j)=n_n+ind;
          else
              elem_mid(2,j)=n_n+ind;
          end
       end
      end % if (for the edge V1-V2)
      
  end % for cyclus over elements
  
  coord_mid=coord_mid(:,1:ind);

end