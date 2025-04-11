function [K, B, WEIGHT] = elastic_stiffness_matrix_3D(ELEM, COORD, shear, bulk, DHatP1, DHatP2, DHatP3, WF)
%--------------------------------------------------------------------------
% elastic_stiffness_matrix_3D assembles the elastic stiffness matrix for a 
% 3D elastoplastic body.
%
% This function computes the elastic stiffness matrix K, the strain-displacement
% matrix B, and the integration weights for each quadrature point. The procedure 
% involves computing the Jacobian of the transformation from the reference 
% element to the current element, evaluating the derivatives of the basis 
% functions in the global coordinate system, assembling the strain-displacement 
% matrix, and finally assembling the elastic stiffness matrix.
%
% INPUT ARGUMENTS:
%   ELEM    - Element connectivity matrix, size: (n_p, n_e)
%             n_p: number of vertices per element, n_e: number of elements.
%   COORD   - Nodal coordinates, size: (3, n_n), where n_n is the number of nodes.
%   shear   - Shear moduli at integration points, size: (1, n_int)
%   bulk    - Bulk moduli at integration points, size: (1, n_int)
%   DHatP1  - Derivatives of local basis functions in the xi_1 direction,
%             size: (n_p, n_q)
%   DHatP2  - Derivatives of local basis functions in the xi_2 direction,
%             size: (n_p, n_q)
%   DHatP3  - Derivatives of local basis functions in the xi_3 direction,
%             size: (n_p, n_q)
%   WF      - Weight factors on the reference element, size: (1, n_q)
%
% OUTPUT:
%   K       - Elastic stiffness matrix, size: (3*n_n, 3*n_n)
%   B       - Strain-displacement matrix, size: (6*n_int, 3*n_n)
%   WEIGHT  - Weight coefficients for each integration point, size: (1, n_int)
%
%--------------------------------------------------------------------------
%
% Auxiliary notation
%
n_n = size(COORD, 2);   % Number of nodes (including midpoints)
n_e = size(ELEM, 2);    % Number of elements
n_p = size(ELEM, 1);    % Number of vertices per element
n_q = length(WF);       % Number of quadrature points per element
n_int = n_e * n_q;      % Total number of integration points

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

%
% Assemble the strain-displacement matrix B.
% B is constructed in a sparse format with size: (6*n_int, 3*n_n)
%

% Build the matrix vB containing values for B.
% n_b is the total number of blocks; each element contributes 18*n_p entries.
n_b = 18 * n_p;
vB = zeros(n_b, n_int);  % Size: (18*n_p, n_int)

% Assign the derivatives into vB corresponding to the proper strain components.
vB(1:18:n_b-17, :) = DPhi1;
vB(10:18:n_b-8, :) = DPhi1;
vB(18:18:n_b, :)   = DPhi1;
vB(4:18:n_b-14, :) = DPhi2;
vB(8:18:n_b-10, :) = DPhi2;
vB(17:18:n_b-1, :) = DPhi2;
vB(6:18:n_b-12, :) = DPhi3;
vB(11:18:n_b-7, :) = DPhi3;
vB(15:18:n_b-3, :) = DPhi3;

% Define indices for assembling the sparse matrix B.
% iB: Row indices for B.
AUX = reshape(1:6*n_int, 6, n_int);
iB = repmat(AUX, 3*n_p, 1);

% jB: Column indices for B.
AUX1 = repmat((1:n_p), 3, 1);    % Size: (3, n_p)
AUX2 = repmat([2; 1; 0], 1, n_p);  % Size: (3, n_p)
% Compute global node indices for each element and degree of freedom.
AUX3 = 3 * ELEM(AUX1(:)', :) - kron(ones(1, n_e), AUX2(:));
jB = kron(AUX3, ones(6, n_q));

% Assemble the sparse strain-displacement matrix B.
B = sparse(iB(:), jB(:), vB(:), 6*n_int, 3*n_n);

%
% Assemble the elastic stress-strain matrix D.
% D is constructed for each integration point and assembled into a sparse matrix.
%

% Define the elastic tensor components.
IOTA = [1; 1; 1; 0; 0; 0];
VOL = IOTA * IOTA'; 
DEV = diag([1, 1, 1, 1/2, 1/2, 1/2]) - VOL / 3;
% ELAST is arranged column-wise for each integration point, size: (36, n_int)
ELAST = 2 * DEV(:) * shear + VOL(:) * bulk;

% Compute the weight coefficients for each integration point.
% Multiply the absolute value of the Jacobian determinant by the quadrature weights.
WEIGHT = abs(DET) .* repmat(WF, 1, n_e);

% Assemble the sparse elastic tensor D.
iD = repmat(AUX, 6, 1);
jD = kron(AUX, ones(6, 1));
vD = ELAST .* (ones(36, 1) * WEIGHT);
D = sparse(iD, jD, vD);

%
% Compute the elastic stiffness matrix.
% K = B' * D * B, size: (3*n_n, 3*n_n)
%
K = B' * D * B;
 
end  % end of function
