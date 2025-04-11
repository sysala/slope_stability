function [Xi, WF] = quadrature_volume_3D(elem_type)
%--------------------------------------------------------------------------
% quadrature_volume_3D defines numerical quadrature rules for volume integration
% in 3D, based on the type of finite element provided.
%
% Depending on the chosen finite element type, this function returns the local 
% coordinates of the quadrature points (Xi) and the corresponding weight factors (WF).
%
% INPUT ARGUMENT:
%   elem_type - string specifying the finite element type. Supported types:
%               'P1' : Linear tetrahedral element (1-point rule)
%               'P2' : Quadratic tetrahedral element (11-point rule)
%               'Q1' : Linear hexahedral element (2x2x2-point rule)
%               'Q2' : Quadratic hexahedral element (3x3x3-point rule)
%
% OUTPUT:
%   Xi       - Local coordinates of quadrature points, size: (3, n_q)
%   WF       - Weight factors corresponding to each quadrature point, size: (1, n_q)
%
%--------------------------------------------------------------------------

switch(elem_type)
    
    case 'P1'
        %------------------------------------------------------------------
        % For the reference tetrahedron with vertices:
        %   [0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]
        % Use a 1-point quadrature rule (n_q = 1)
        %------------------------------------------------------------------
        Xi = [1/4; 1/4; 1/4];
        WF = 1/6;
        
    case 'P2'
        %------------------------------------------------------------------
        % For the reference tetrahedron with vertices:
        %   [0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]
        % Use an 11-point quadrature rule (n_q = 11)
        %------------------------------------------------------------------
        Xi = [1/4, 0.0714285714285714, 0.785714285714286, 0.0714285714285714, 0.0714285714285714, 0.399403576166799, 0.100596423833201, 0.100596423833201, 0.399403576166799, 0.399403576166799, 0.100596423833201;
              1/4, 0.0714285714285714, 0.0714285714285714, 0.785714285714286, 0.0714285714285714, 0.100596423833201, 0.399403576166799, 0.100596423833201, 0.399403576166799, 0.100596423833201, 0.399403576166799;
              1/4, 0.0714285714285714, 0.0714285714285714, 0.0714285714285714, 0.785714285714286, 0.100596423833201, 0.100596423833201, 0.399403576166799, 0.100596423833201, 0.399403576166799, 0.399403576166799];
        WF = [-0.013155555555555, 0.007622222222222*ones(1,4), 0.024888888888888*ones(1,6)];
        
    case 'Q1'
        %------------------------------------------------------------------
        % For the reference cube with vertices:
        %   [-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
        %   [-1, -1, 1],  [1, -1, 1],  [1, 1, 1],  [-1, 1, 1]
        % Use a (2x2x2)-point quadrature rule (n_q = 8)
        %------------------------------------------------------------------
        pt = 1/sqrt(3);
        Xi = [ -pt, -pt, -pt, -pt,  pt,  pt,  pt,  pt;
               -pt, -pt,  pt,  pt, -pt, -pt,  pt,  pt;
               -pt,  pt, -pt,  pt, -pt,  pt, -pt,  pt ];
        WF = [1, 1, 1, 1, 1, 1, 1, 1];
        
    case 'Q2'
        %------------------------------------------------------------------
        % For the reference cube with vertices:
        %   [-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
        %   [-1, -1, 1],  [1, -1, 1],  [1, 1, 1],  [-1, 1, 1]
        % Use a (3x3x3)-point quadrature rule (n_q = 27)
        %------------------------------------------------------------------
        pt = sqrt(3/5);
        Xi = [ -pt, -pt, -pt, -pt, -pt, -pt, -pt, -pt, -pt,   0,   0,   0,   0,   0,   0,   0,   0,   0,  pt,  pt,  pt,  pt,  pt,  pt,  pt,  pt,  pt;
               -pt, -pt, -pt,   0,   0,   0,  pt,  pt,  pt, -pt, -pt, -pt,   0,   0,   0,  pt,  pt,  pt, -pt, -pt, -pt,   0,   0,   0,  pt,  pt,  pt;
               -pt,   0,  pt, -pt,   0,  pt, -pt,   0,  pt, -pt,   0,  pt, -pt,   0,  pt, -pt,   0,  pt, -pt,   0,  pt, -pt,   0,  pt, -pt,   0,  pt ];
        % The weights WF are computed as a weighted sum of several terms.
        WF = ((5/9)^3)* [1 0 1 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 1 0 1] + ...
             ((5/9)^2)*(8/9)* [0 1 0 1 0 1 0 1 0 1 0 1 0 0 0 1 0 1 0 1 0 1 0 1 0 1 0] + ...     
             ((8/9)^2)*(5/9)* [0 0 0 0 1 0 0 0 0 0 1 0 1 0 1 0 1 0 0 0 0 0 1 0 0 0 0] + ...
             ((8/9)^3)* [0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0];
        
    otherwise
        disp('Bad choice of element type');
end
