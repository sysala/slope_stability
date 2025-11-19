function [coord_mid,elem_mid,surf]=midpoints_P4(coord,elem)

%  The aim of this function is to complete basic mesh data (coord, elem)
%  with infomation about midpoints and surface nodes. P4 elements are
%  considered.
% 
%  input data:
%    coord - coordinates of the vertices, size(coord)=(2,n_n) where n_n is 
%            a number of vertices
%    elem - 3 x n_e array containing numbers of nodes defining each
%           element, n_e = number of elements. We keep anticlockwise
%           ordering of the nodes creating a particular element. 
%
%  output data:
%    coord_mid - coordinates of midpoints
%    elem_mid - 10 x n_e array containing numbers of midpoints in each
%               element
%    surf    -  5 x n_surf arrays containg numbers of nodal 
%               points defining edges of the domain boundary 
%
% ======================================================================
%
  
  % numbers of elements and vertices
  n_e=size(elem,2);
  n_n=size(coord,2);
      
  % predefinition of unknown arrays
  coord_mid=zeros(2,10*n_e);
  elem_mid=zeros(12,n_e);
  surf=zeros(5,n_e);
  
  % for cyclus over elements
  ind=0;   % enlarging index specifying midpoints
  ind_s=0; % enlarging index specifying surf
  for i=1:n_e
      
      % vertices defining the i-th element 
      V1=elem(1,i);
      V2=elem(2,i);
      V3=elem(3,i);
      
      % creation of the midpoints which do not belong on edges
      coord_mid(:,ind+1)=coord(:,V1)/2+coord(:,V2)/4+coord(:,V3)/4; 
      elem_mid(10,i)=n_n+ind+1;   
      %
      coord_mid(:,ind+2)=coord(:,V1)/4+coord(:,V2)/2+coord(:,V3)/4; 
      elem_mid(11,i)=n_n+ind+2;  
      %
      coord_mid(:,ind+3)=coord(:,V1)/4+coord(:,V2)/4+coord(:,V3)/2; 
      elem_mid(12,i)=n_n+ind+3;   
      %
      ind=ind+3;
      
      % analysis of the edge V1-V2
      if elem_mid(1,i)==0
       % creation of new midpoints lying on the edge V1-V2 
       coord_mid(:,ind+1)=(coord(:,V1)+coord(:,V2))/2; 
       elem_mid(1,i)=n_n+ind+1;   
       %
       coord_mid(:,ind+2)=3*coord(:,V1)/4+coord(:,V2)/4; 
       elem_mid(4,i)=n_n+ind+2;   
       %
       coord_mid(:,ind+3)=coord(:,V1)/4+3*coord(:,V2)/4; 
       elem_mid(5,i)=n_n+ind+3;   
       % finding the adjacent element j to i which contains the edge V1-V2
       [row1,col1]=find(elem==V1);
       [row2,col2]=find(elem==V2);
       j=setdiff(intersect(col1,col2),i);
       if ~isempty(j)
          % This case means that the edge V1-V2 is the intersection of the
          % elements i and j.
          v=elem(:,j);
          if V2==v(1)
              elem_mid(1,j)=n_n+ind+1;
              elem_mid(4,j)=n_n+ind+3;
              elem_mid(5,j)=n_n+ind+2;
          elseif V2==v(2)
              elem_mid(2,j)=n_n+ind+1;
              elem_mid(6,j)=n_n+ind+3;
              elem_mid(7,j)=n_n+ind+2;
          else
              elem_mid(3,j)=n_n+ind+1;
              elem_mid(8,j)=n_n+ind+3;
              elem_mid(9,j)=n_n+ind+2;
          end
       else
          ind_s=ind_s+1;  
          surf(:,ind_s)=[V1; V2; n_n+ind+1; n_n+ind+2; n_n+ind+3];
       end
       ind=ind+3;
      end % if (for the edge V1-V2)
      
     % analysis of the edge V2-V3
      if elem_mid(2,i)==0
       % creation of new midpoints lying on the edge V2-V3 
       coord_mid(:,ind+1)=(coord(:,V2)+coord(:,V3))/2; 
       elem_mid(2,i)=n_n+ind+1;   
       %
       coord_mid(:,ind+2)=3*coord(:,V2)/4+coord(:,V3)/4; 
       elem_mid(6,i)=n_n+ind+2;   
       %
       coord_mid(:,ind+3)=coord(:,V2)/4+3*coord(:,V3)/4; 
       elem_mid(7,i)=n_n+ind+3;   
       % finding the adjacent element j to i which contains the edge V2-V3
       [row1,col1]=find(elem==V2);
       [row2,col2]=find(elem==V3);
       j=setdiff(intersect(col1,col2),i);
       if ~isempty(j)
          % This case means that the edge V2-V3 is the intersection of the
          % elements i and j.
          v=elem(:,j);
          if V3==v(1)
              elem_mid(1,j)=n_n+ind+1;
              elem_mid(4,j)=n_n+ind+3;
              elem_mid(5,j)=n_n+ind+2;
          elseif V3==v(2)
              elem_mid(2,j)=n_n+ind+1;
              elem_mid(6,j)=n_n+ind+3;
              elem_mid(7,j)=n_n+ind+2;
          else
              elem_mid(3,j)=n_n+ind+1;
              elem_mid(8,j)=n_n+ind+3;
              elem_mid(9,j)=n_n+ind+2;
          end
       else
          ind_s=ind_s+1;  
          surf(:,ind_s)=[V1; V2; n_n+ind+1; n_n+ind+2; n_n+ind+3];
       end
       ind=ind+3;
      end % if (for the edge V2-V3)
      
      % analysis of the edge V3-V1
      if elem_mid(3,i)==0
       % creation of new midpoints lying on the edge V3-V1 
       coord_mid(:,ind+1)=(coord(:,V3)+coord(:,V1))/2; 
       elem_mid(3,i)=n_n+ind+1;   
       %
       coord_mid(:,ind+2)=3*coord(:,V3)/4+coord(:,V1)/4; 
       elem_mid(8,i)=n_n+ind+2;   
       %
       coord_mid(:,ind+3)=coord(:,V3)/4+3*coord(:,V1)/4; 
       elem_mid(9,i)=n_n+ind+3;   
       % finding the adjacent element j to i which contains the edge V3-V1
       [row1,col1]=find(elem==V3);
       [row2,col2]=find(elem==V1);
       j=setdiff(intersect(col1,col2),i);
       if ~isempty(j)
          % This case means that the edge V3-V1 is the intersection of the
          % elements i and j.
          v=elem(:,j);
          if V1==v(1)
              elem_mid(1,j)=n_n+ind+1;
              elem_mid(4,j)=n_n+ind+3;
              elem_mid(5,j)=n_n+ind+2;
          elseif V1==v(2)
              elem_mid(2,j)=n_n+ind+1;
              elem_mid(6,j)=n_n+ind+3;
              elem_mid(7,j)=n_n+ind+2;
          else
              elem_mid(3,j)=n_n+ind+1;
              elem_mid(8,j)=n_n+ind+3;
              elem_mid(9,j)=n_n+ind+2;
          end
       else
          ind_s=ind_s+1;  
          surf(:,ind_s)=[V1; V2; n_n+ind+1; n_n+ind+2; n_n+ind+3];
       end
       ind=ind+3;
      end % if (for the edge V3-V1)
      
  end % for cyclus over elements
  
  coord_mid=coord_mid(:,1:ind);
  surf=surf(:,1:ind_s);
 
end