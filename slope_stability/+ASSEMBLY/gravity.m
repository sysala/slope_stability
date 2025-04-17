function gamma=gravity(gamma_sat,gamma_unsat,coord,elem,HatP)

% =========================================================================
%
%  This function determines the specific weight at integration points 
%  depending on a given saturation curve. This curve is a part of this
%  code.
%
%  input data:
%    gamma_sat -  1 x n_int array of specific weights at integration points
%                 for saturated materials
%    gamma_unsat- 1 x n_int array of specific weights at integration points
%                 for unsaturated materials
%    coord - coordinates of the nodes, size(coord)=(2,n_n) where n_n is a
%            number of nodes
%    elem - n_p x n_e array containing numbers of vertices defining each
%           element, n_e = number of elements
%    HatP - values of basis functions at the quadrature points,
%           size(HatP)=(n_p,n_q)
% ======================================================================
%

%
% Auxilliary notation
%

  n_e=size(elem,2);     % number of elements
  n_p=size(elem,1);     % number of vertices per element
  n_q=size(HatP,2);     % number of quadrature points
  n_int = n_e*n_q ;     % total number of integrations points
     
%
% Coordinates of integration points 
% 

  % extension of the input arrays HatP by replication
  % size(HatPhi)=(n_p, n_int)
  HatPhi=repmat(HatP,1,n_e);
  
  % coordinates of nodes defining each element
  % size(COORDe1)=size(COORDe2)=(n_p, n_e)
  COORDe1=reshape(coord(1,elem(:)),n_p,n_e);
  COORDe2=reshape(coord(2,elem(:)),n_p,n_e);
  
  % coordinates of nodes around each integration point
  % size(COORDint1)=size(COORDint2)=(n_p, n_int)
  COORDint1=kron(COORDe1,ones(1,n_q)); 
  COORDint2=kron(COORDe2,ones(1,n_q)); 
  
  % coordinates of integration points: 
  %  size(C_int1)=size(C_int2)=(1, n_int)
  C_int1=sum(COORDint1.*HatPhi); C_int2=sum(COORDint2.*HatPhi); 

%  
% given saturation curve (level of porous water flow)
%
  function y=f(x)
     if x<=44
         y=59-(59-55)*(x-0)/(44-0);
     elseif x<=116
         y=55-(55-39)*(x-44)/(116-44);
     elseif x<=149
         y=39-(39-32)*(x-116)/(149-116);         
     elseif x<=165
         y=32-(32-27)*(x-149)/(165-149);  
     elseif x<=194
         y=27-(27-24)*(x-165)/(194-165);           
     elseif x<=232
         y=24-(24-20)*(x-194)/(232-194);    
     else
         y=20;
     end
  end

%
% Assembling of the output array gamma
% 
  gamma=zeros(1,n_int);
  for i=1:n_int
      x=C_int1(i); y=C_int2(i);
      if y<=f(x)
          gamma(i)=gamma_sat(i);
      else
          gamma(i)=gamma_unsat(i);
      end
  end       
  
end % function