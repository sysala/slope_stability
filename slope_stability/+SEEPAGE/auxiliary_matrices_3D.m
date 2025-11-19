
function [B,C,WEIGHT]=auxiliary_matrices_3D...
                                  (ELEM,COORD,HatP,DHatP1,DHatP2,DHatP3,WF)
 
% =========================================================================
%
%  This function assembles auxiliary matrices for the construction of
%  linearized stiffness matrices and rhs vectors within Newton-like solvers
%
%  input data:
%    ELEM - n_p x n_e array containing numbers of nodes defining each
%           element, n_e = number of elements
%    COORD - coordinates of the nodes, size(coord)=(2,n_n) where n_n is a
%            number of nodes
%    HatP   - values of basis functions at the quadrature points, 
%             size(HatP)=(n_p,n_q)
%    DHatP1  - Derivatives of local basis functions in the xi_1 direction,
%             size: (n_p, n_q)
%    DHatP2  - Derivatives of local basis functions in the xi_2 direction,
%             size: (n_p, n_q)
%    DHatP3  - Derivatives of local basis functions in the xi_3 direction,
%             size: (n_p, n_q)
%        WF - weight factors on the reference element, size(WF)=(1,n_q)
%
%  output data:
%    B - a matrix representing the gradient operator, size(B)=(3*n_int,n_n)
%    C - a matrix representing function values, size(C)=(n_int,n_n)
%    WEIGHT - the weight coefficients for each quadrature point, 
%             size(WEIGHT)=(1, n_int)
%
% ======================================================================
%

%
% Auxilliary notation
%
n_n=size(COORD,2);    % number of nodes including midpoints
n_e=size(ELEM,2);     % number of elements
n_p=size(ELEM,1);     % number of vertices per element
n_q=length(WF);       % number of quadrature points
n_int = n_e*n_q ;     % total number of integrations points
     
%
% Compute the Jacobian, its determinant, inverse, and the derivatives of 
% local basis functions with respect to the global coordinates.
%

% Replicate local basis function derivatives to cover all integration points.
% Resulting size: (n_p, n_int)
DHatPhi1 = repmat(DHatP1, 1, n_e);
DHatPhi2 = repmat(DHatP2, 1, n_e);
DHatPhi3 = repmat(DHatP3, 1, n_e);

% Extract nodal coordinates for each element.
% Each COORDe* has size: (n_p, n_e)
COORDe1 = reshape(COORD(1, ELEM(:)), n_p, n_e);
COORDe2 = reshape(COORD(2, ELEM(:)), n_p, n_e);
COORDe3 = reshape(COORD(3, ELEM(:)), n_p, n_e);

% Map nodal coordinates to each integration point.
% Each COORDint* has size: (n_p, n_int)
COORDint1 = kron(COORDe1, ones(1, n_q));
COORDint2 = kron(COORDe2, ones(1, n_q));
COORDint3 = kron(COORDe3, ones(1, n_q));

% Compute components of the Jacobian at each integration point.
J11 = sum(COORDint1 .* DHatPhi1);
J12 = sum(COORDint2 .* DHatPhi1);
J13 = sum(COORDint3 .* DHatPhi1);
J21 = sum(COORDint1 .* DHatPhi2);
J22 = sum(COORDint2 .* DHatPhi2);
J23 = sum(COORDint3 .* DHatPhi2);
J31 = sum(COORDint1 .* DHatPhi3);
J32 = sum(COORDint2 .* DHatPhi3);
J33 = sum(COORDint3 .* DHatPhi3);

% Compute determinant of the Jacobian.
DET = J11 .* (J22 .* J33 - J32 .* J23) - J12 .* (J21 .* J33 - J23 .* J31) + J13 .* (J21 .* J32 - J22 .* J31);

% Compute components of the inverse Jacobian.
Jinv11 =  (J22 .* J33 - J23 .* J32) ./ DET; 
Jinv12 = -(J12 .* J33 - J13 .* J32) ./ DET; 
Jinv13 =  (J12 .* J23 - J13 .* J22) ./ DET; 

Jinv21 = -(J21 .* J33 - J23 .* J31) ./ DET; 
Jinv22 =  (J11 .* J33 - J13 .* J31) ./ DET; 
Jinv23 = -(J11 .* J23 - J13 .* J21) ./ DET; 

Jinv31 =  (J21 .* J32 - J22 .* J31) ./ DET; 
Jinv32 = -(J11 .* J32 - J12 .* J31) ./ DET; 
Jinv33 =  (J11 .* J22 - J12 .* J21) ./ DET; 

% Transform derivatives of basis functions to global coordinates.
% Each DPhi* has size: (n_p, n_int)
DPhi1 = repmat(Jinv11, n_p, 1) .* DHatPhi1 + repmat(Jinv12, n_p, 1) .* DHatPhi2 + repmat(Jinv13, n_p, 1) .* DHatPhi3;
DPhi2 = repmat(Jinv21, n_p, 1) .* DHatPhi1 + repmat(Jinv22, n_p, 1) .* DHatPhi2 + repmat(Jinv23, n_p, 1) .* DHatPhi3;
DPhi3 = repmat(Jinv31, n_p, 1) .* DHatPhi1 + repmat(Jinv32, n_p, 1) .* DHatPhi2 + repmat(Jinv33, n_p, 1) .* DHatPhi3;
  
% Compute the weight coefficients for each integration point.
% Multiply the absolute value of the Jacobian determinant by the quadrature weights.
WEIGHT = abs(DET) .* repmat(WF, 1, n_e);
  
%
% assembling of the gradient-pressure matrix B
% size(B)=(3*n_int, n_n)
%   
  
  % values of the gradient-pressure matrix B
  vB=zeros(3*n_p,n_int);        % size(vB)=(3*n_p,n_int)
  vB(1:3:3*n_p-2,:)=DPhi1; 
  vB(2:3:3*n_p-1,:)=DPhi2; 
  vB(3:3:3*n_p  ,:)=DPhi3; 

  % i-th and j-th indices of B: size(iB)=size(jB)=(3*n_p,n_int)
  AUX=reshape(1:3*n_int,3,n_int);
  iB=repmat(AUX,n_p,1);    
  jB=kron(ELEM,ones(3,n_q));
  
  % the sparse gradient-pressure matrix B
  B=sparse(iB(:),jB(:),vB(:), 3*n_int,n_n);  
  
%
% assembling of the matrix C
% size(C)=(n_int, n_n)
%   
  
  % values of the matrix C
  vC=repmat(HatP,1,n_e);        % size(vC)=(n_p,n_int)
 
  % i-th and j-th indices of C: size(iC)=size(jC)=(n_p,n_int)
  iC=repmat(1:n_int,n_p,1);    
  jC=kron(ELEM,ones(1,n_q));
  
  % the sparse matrix C
  C=sparse(iC(:),jC(:),vC(:), n_int,n_n);    
 
end  % end of function