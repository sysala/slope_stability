function surf= surf_P1_2D(elem)

%  The aim of this function is to complete the array elem
%  with numbers of surface nodes. P1 elements are considered.
%
%  input data:
%    elem - 3 x n_e array containing numbers of nodes defining each
%           element, n_e = number of elements. We keep anticlockwise
%           ordering of the nodes creating a particular element. 
%
%  output data:
%    surf    -  2 x n_surf arrays containg numbers of nodal 
%               points defining edges of the domain boundary 
%
% ======================================================================
%
  
  % numbers of elements and vertices
  n_e=size(elem,2);
      
  % predefinition of unknown arrays
  surf=zeros(2,n_e);
  
  % for cyclus over elements
  ind_s=0; % enlarging index specifying surf
  for i=1:n_e
      
      % vertices defining the i-th element 
      V1=elem(1,i);
      V2=elem(2,i);
      V3=elem(3,i);
      
      % analysis of the edge V2-V3
      % finding the adjacent element j to i which contains the edge V2-V3
       [row1,col1]=find(elem==V2);
       [row2,col2]=find(elem==V3);
       j=setdiff(intersect(col1,col2),i);
       if isempty(j)
          ind_s=ind_s+1;  
          surf(:,ind_s)=[V3; V2];
       end
       
      % analysis of the edge V3-V1           
      % finding the adjacent element j to i which contains the edge V3-V1
       [row1,col1]=find(elem==V3);
       [row2,col2]=find(elem==V1);
       j=setdiff(intersect(col1,col2),i);
       if isempty(j)
          ind_s=ind_s+1;  
          surf(:,ind_s)=[V3; V2];
       end
           
      % analysis of the edge V1-V2
      % finding the adjacent element j to i which contains the edge V1-V2
       [row1,col1]=find(elem==V1);
       [row2,col2]=find(elem==V2);
       j=setdiff(intersect(col1,col2),i);
       if isempty(j)
          ind_s=ind_s+1;  
          surf(:,ind_s)=[V3; V2];
       end
      
  end % for cyclus over elements
  
  surf=surf(:,1:ind_s);

end