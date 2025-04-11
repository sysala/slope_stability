function [HatP, DHatP1, DHatP2, DHatP3] = local_basis_volume_3D(elem_type, Xi)
%--------------------------------------------------------------------------
% local_basis_volume_3D evaluates local basis functions and their derivatives
% for volume integration in 3D, depending on the chosen finite element type.
%
% This function computes the values of the basis functions (HatP) and their
% derivatives (DHatP1, DHatP2, DHatP3) at prescribed quadrature points given
% by the matrix Xi. The evaluation is performed according to the type of
% Lagrange finite elements specified by elem_type.
%
% INPUT ARGUMENTS:
%   elem_type - string specifying the finite element type. Supported types:
%               'P1' : Linear tetrahedral element.
%               'P2' : Quadratic tetrahedral element.
%               'Q1' : Linear hexahedral element.
%               'Q2' : Quadratic hexahedral element.
%   Xi        - Local coordinates of the quadrature points, size: (3, n_q),
%               where n_q is the number of quadrature points.
%
% OUTPUT:
%   HatP    - Values of the basis functions at the quadrature points,
%             size: (n_p, n_q), where n_p is the number of basis functions.
%   DHatP1  - Derivatives of the basis functions in the xi_1 direction,
%             size: (n_p, n_q)
%   DHatP2  - Derivatives of the basis functions in the xi_2 direction,
%             size: (n_p, n_q)
%   DHatP3  - Derivatives of the basis functions in the xi_3 direction,
%             size: (n_p, n_q)
%
%--------------------------------------------------------------------------

% Extract individual coordinates from Xi for clarity.
xi_1 = Xi(1, :); 
xi_2 = Xi(2, :); 
xi_3 = Xi(3, :);

switch(elem_type)
    case 'P1'
        %------------------------------------------------------------------
        % For the reference tetrahedron with vertices:
        %   [0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]
        % n_p = 4, n_q = length(xi_1)
        %------------------------------------------------------------------
        % Basis functions defined using barycentric coordinates.
        HatP = [1 - xi_1 - xi_2 - xi_3; xi_1; xi_2; xi_3];
        % Derivatives of the basis functions are constant.
        DHatP1 = [-1; 1; 0; 0];
        DHatP2 = [-1; 0; 1; 0];
        DHatP3 = [-1; 0; 0; 1];
        
    case 'P2'
        %------------------------------------------------------------------
        % For the reference tetrahedron with vertices:
        %   [0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]
        % and mid-edge nodes:
        %   [1/2, 0, 0], [1/2, 1/2, 0], [0, 1/2, 0],
        %   [1/2, 0, 1/2], [0, 1/2, 1/2], [0, 0, 1/2]
        % n_p = 10, n_q = length(xi_1)
        %------------------------------------------------------------------
        % Barycentric coordinate for the first vertex.
        xi_0 = 1 - xi_1 - xi_2 - xi_3;
        n_q = length(xi_1);
        % Quadratic basis functions in barycentric coordinates.
        HatP = [ xi_0 .* (2*xi_0 - 1);
                 xi_1 .* (2*xi_1 - 1);
                 xi_2 .* (2*xi_2 - 1);
                 xi_3 .* (2*xi_3 - 1);
                 4 * xi_0 .* xi_1;
                 4 * xi_1 .* xi_2;
                 4 * xi_0 .* xi_2;
                 4 * xi_1 .* xi_3;
                 4 * xi_2 .* xi_3;
                 4 * xi_0 .* xi_3 ];
        % Derivatives with respect to xi_1.
        DHatP1 = [ -4*xi_0 + 1;
                    4*xi_1 - 1;
                    zeros(1, n_q);
                    zeros(1, n_q);
                    4*(xi_0 - xi_1);
                    4*xi_2;
                    -4*xi_2;
                    4*xi_3;
                    zeros(1, n_q);
                    -4*xi_3 ];
        % Derivatives with respect to xi_2.
        DHatP2 = [ -4*xi_0 + 1;
                    zeros(1, n_q);
                    4*xi_2 - 1;
                    zeros(1, n_q);
                    -4*xi_1;
                    4*xi_1;
                    4*(xi_0 - xi_2);
                    zeros(1, n_q);
                    4*xi_3;
                    -4*xi_3 ];
        % Derivatives with respect to xi_3.
        DHatP3 = [ -4*xi_0 + 1;
                    zeros(1, n_q);
                    zeros(1, n_q);
                    4*xi_3 - 1;
                    -4*xi_1;
                    zeros(1, n_q);
                    -4*xi_2;
                    4*xi_1;
                    4*xi_2;
                    4*(xi_0 - xi_3) ];
        
    case 'Q1'
        %------------------------------------------------------------------
        % For the reference cube with vertices:
        %   [-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
        %   [-1, -1, 1],  [1, -1, 1],  [1, 1, 1],  [-1, 1, 1]
        % n_p = 8, n_q = length(xi_1)
        %------------------------------------------------------------------
        % Trilinear basis functions for the hexahedral element.
        HatP = [ (1 - xi_1).*(1 - xi_2).*(1 - xi_3) / 8;
                 (1 + xi_1).*(1 - xi_2).*(1 - xi_3) / 8;
                 (1 + xi_1).*(1 + xi_2).*(1 - xi_3) / 8;
                 (1 - xi_1).*(1 + xi_2).*(1 - xi_3) / 8;
                 (1 - xi_1).*(1 - xi_2).*(1 + xi_3) / 8;
                 (1 + xi_1).*(1 - xi_2).*(1 + xi_3) / 8;
                 (1 + xi_1).*(1 + xi_2).*(1 + xi_3) / 8;
                 (1 - xi_1).*(1 + xi_2).*(1 + xi_3) / 8 ];
        % Derivatives in the xi_1 direction.
        DHatP1 = [ - (1 - xi_2).*(1 - xi_3) / 8;
                    (1 - xi_2).*(1 - xi_3) / 8;
                    (1 + xi_2).*(1 - xi_3) / 8;
                   - (1 + xi_2).*(1 - xi_3) / 8;
                   - (1 - xi_2).*(1 + xi_3) / 8;
                    (1 - xi_2).*(1 + xi_3) / 8;
                    (1 + xi_2).*(1 + xi_3) / 8;
                   - (1 + xi_2).*(1 + xi_3) / 8 ];
        % Derivatives in the xi_2 direction.
        DHatP2 = [ - (1 - xi_1).*(1 - xi_3) / 8;
                   - (1 + xi_1).*(1 - xi_3) / 8;
                    (1 + xi_1).*(1 - xi_3) / 8;
                    (1 - xi_1).*(1 - xi_3) / 8;
                   - (1 - xi_1).*(1 + xi_3) / 8;
                   - (1 + xi_1).*(1 + xi_3) / 8;
                    (1 + xi_1).*(1 + xi_3) / 8;
                    (1 - xi_1).*(1 + xi_3) / 8 ];
        % Derivatives in the xi_3 direction.
        DHatP3 = [ - (1 - xi_1).*(1 - xi_2) / 8;
                   - (1 + xi_1).*(1 - xi_2) / 8;
                   - (1 + xi_1).*(1 + xi_2) / 8;
                   - (1 - xi_1).*(1 + xi_2) / 8;
                    (1 - xi_1).*(1 - xi_2) / 8;
                    (1 + xi_1).*(1 - xi_2) / 8;
                    (1 + xi_1).*(1 + xi_2) / 8;
                    (1 - xi_1).*(1 + xi_2) / 8 ];
        
    case 'Q2'
        %------------------------------------------------------------------
        % For the reference cube with vertices:
        %   [-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
        %   [-1, -1, 1],  [1, -1, 1],  [1, 1, 1],  [-1, 1, 1]
        % and additional mid-side nodes,
        % n_p = 20, n_q = length(xi_1)
        %------------------------------------------------------------------
        % Quadratic basis functions for the hexahedral element.
        HatP = [ (1 - xi_1).*(1 - xi_2).*(1 - xi_3).*(-2 - xi_1 - xi_2 - xi_3) / 8;
                 (1 + xi_1).*(1 - xi_2).*(1 - xi_3).*(-2 + xi_1 - xi_2 - xi_3) / 8;
                 (1 + xi_1).*(1 + xi_2).*(1 - xi_3).*(-2 + xi_1 + xi_2 - xi_3) / 8;
                 (1 - xi_1).*(1 + xi_2).*(1 - xi_3).*(-2 - xi_1 + xi_2 - xi_3) / 8;
                 (1 - xi_1).*(1 - xi_2).*(1 + xi_3).*(-2 - xi_1 - xi_2 + xi_3) / 8;
                 (1 + xi_1).*(1 - xi_2).*(1 + xi_3).*(-2 + xi_1 - xi_2 + xi_3) / 8;
                 (1 + xi_1).*(1 + xi_2).*(1 + xi_3).*(-2 + xi_1 + xi_2 + xi_3) / 8;
                 (1 - xi_1).*(1 + xi_2).*(1 + xi_3).*(-2 - xi_1 + xi_2 + xi_3) / 8;
                 (1 - xi_1.^2).*(1 - xi_2).*(1 - xi_3) / 4;
                 (1 + xi_1).*(1 - xi_2.^2).*(1 - xi_3) / 4;
                 (1 - xi_1.^2).*(1 + xi_2).*(1 - xi_3) / 4;
                 (1 - xi_1).*(1 - xi_2.^2).*(1 - xi_3) / 4;
                 (1 - xi_1.^2).*(1 - xi_2).*(1 + xi_3) / 4;
                 (1 + xi_1).*(1 - xi_2.^2).*(1 + xi_3) / 4;
                 (1 - xi_1.^2).*(1 + xi_2).*(1 + xi_3) / 4;
                 (1 - xi_1).*(1 - xi_2.^2).*(1 + xi_3) / 4;
                 (1 - xi_1).*(1 - xi_2).*(1 - xi_3.^2) / 4;
                 (1 + xi_1).*(1 - xi_2).*(1 - xi_3.^2) / 4;
                 (1 + xi_1).*(1 + xi_2).*(1 - xi_3.^2) / 4;
                 (1 - xi_1).*(1 + xi_2).*(1 - xi_3.^2) / 4 ];
             
        % Derivatives of quadratic basis functions for Q2 elements.
        % The derivatives are computed component-wise.
        DHatP1 = [ (1 - xi_2).*(1 - xi_3).*( 1 + 2*xi_1 + xi_2 + xi_3) / 8;
                   (1 - xi_2).*(1 - xi_3).*(-1 + 2*xi_1 - xi_2 - xi_3) / 8;
                   (1 + xi_2).*(1 - xi_3).*(-1 + 2*xi_1 + xi_2 - xi_3) / 8;
                   (1 + xi_2).*(1 - xi_3).*( 1 + 2*xi_1 - xi_2 + xi_3) / 8;
                   (1 - xi_2).*(1 + xi_3).*( 1 + 2*xi_1 + xi_2 - xi_3) / 8;
                   (1 - xi_2).*(1 + xi_3).*(-1 + 2*xi_1 - xi_2 + xi_3) / 8;
                   (1 + xi_2).*(1 + xi_3).*(-1 + 2*xi_1 + xi_2 + xi_3) / 8;
                   (1 + xi_2).*(1 + xi_3).*( 1 + 2*xi_1 - xi_2 - xi_3) / 8;
                   -xi_1 .* (1 - xi_2).*(1 - xi_3) / 2;
                   (1 - xi_2.^2).*(1 - xi_3) / 4;
                   -xi_1 .* (1 + xi_2).*(1 - xi_3) / 2;
                   -(1 - xi_2.^2).*(1 - xi_3) / 4;
                   -xi_1 .* (1 - xi_2).*(1 + xi_3) / 2;
                   (1 - xi_2.^2).*(1 + xi_3) / 4;
                   -xi_1 .* (1 + xi_2).*(1 + xi_3) / 2;
                   -(1 - xi_2.^2).*(1 + xi_3) / 4;
                   -(1 - xi_2).*(1 - xi_3.^2) / 4;
                   (1 - xi_2).*(1 - xi_3.^2) / 4;
                   (1 + xi_2).*(1 - xi_3.^2) / 4;
                   -(1 + xi_2).*(1 - xi_3.^2) / 4 ];
               
        DHatP2 = [ (1 - xi_1).*(1 - xi_3).*( 1 + xi_1 + 2*xi_2 + xi_3) / 8;
                   (1 + xi_1).*(1 - xi_3).*( 1 - xi_1 + 2*xi_2 + xi_3) / 8;
                   (1 + xi_1).*(1 - xi_3).*(-1 + xi_1 + 2*xi_2 - xi_3) / 8;
                   (1 - xi_1).*(1 - xi_3).*(-1 - xi_1 + 2*xi_2 - xi_3) / 8;
                   (1 - xi_1).*(1 + xi_3).*( 1 + xi_1 + 2*xi_2 - xi_3) / 8;
                   (1 + xi_1).*(1 + xi_3).*( 1 - xi_1 + 2*xi_2 - xi_3) / 8;
                   (1 + xi_1).*(1 + xi_3).*(-1 + xi_1 + 2*xi_2 + xi_3) / 8;
                   (1 - xi_1).*(1 + xi_3).*(-1 - xi_1 + 2*xi_2 + xi_3) / 8;
                   -(1 - xi_1.^2).*(1 - xi_3) / 4;
                   -(1 + xi_1).*xi_2.*(1 - xi_3) / 2;
                   (1 - xi_1.^2).*(1 - xi_3) / 4;
                   -(1 - xi_1).*xi_2.*(1 - xi_3) / 2;
                   -(1 - xi_1.^2).*(1 + xi_3) / 4;
                   -(1 + xi_1).*xi_2.*(1 + xi_3) / 2;
                   (1 - xi_1.^2).*(1 + xi_3) / 4;
                   -(1 - xi_1).*xi_2.*(1 + xi_3) / 2;
                   -(1 - xi_1).*(1 - xi_3.^2) / 4;
                   -(1 + xi_1).*(1 - xi_3.^2) / 4;
                   (1 + xi_1).*(1 - xi_3.^2) / 4;
                   (1 - xi_1).*(1 - xi_3.^2) / 4 ];
               
        DHatP3 = [ (1 - xi_1).*(1 - xi_2).*( 1 + xi_1 + xi_2 + 2*xi_3) / 8;
                   (1 + xi_1).*(1 - xi_2).*( 1 - xi_1 + xi_2 + 2*xi_3) / 8;
                   (1 + xi_1).*(1 + xi_2).*( 1 - xi_1 - xi_2 + 2*xi_3) / 8;
                   (1 - xi_1).*(1 + xi_2).*( 1 + xi_1 - xi_2 + 2*xi_3) / 8;
                   (1 - xi_1).*(1 - xi_2).*(-1 - xi_1 - xi_2 + 2*xi_3) / 8;
                   (1 + xi_1).*(1 - xi_2).*(-1 + xi_1 - xi_2 + 2*xi_3) / 8;
                   (1 + xi_1).*(1 + xi_2).*(-1 + xi_1 + xi_2 + 2*xi_3) / 8;
                   (1 - xi_1).*(1 + xi_2).*(-1 - xi_1 + xi_2 + 2*xi_3) / 8;
                   -(1 - xi_1.^2).*(1 - xi_2) / 4;
                   -(1 + xi_1).*(1 - xi_2.^2) / 4;
                   -(1 - xi_1.^2).*(1 + xi_2) / 4;
                   -(1 - xi_1).*(1 - xi_2.^2) / 4;
                    (1 - xi_1.^2).*(1 - xi_2) / 4;
                    (1 + xi_1).*(1 - xi_2.^2) / 4;
                    (1 - xi_1.^2).*(1 + xi_2) / 4;
                    (1 - xi_1).*(1 - xi_2.^2) / 4;
                   -(1 - xi_1).*(1 - xi_2).*xi_3 / 2;
                   -(1 + xi_1).*(1 - xi_2).*xi_3 / 2;
                   -(1 + xi_1).*(1 + xi_2).*xi_3 / 2;
                   -(1 - xi_1).*(1 + xi_2).*xi_3 / 2 ];
        
    otherwise
        disp('Bad choice of element type');
end
