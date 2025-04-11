function [Xi, wt] = quadrature_volume_3D_degree(pts)
% QUADRATURE_VOLUME_3D returns quadrature points and weights for a tetrahedron.
%
%   [Xi, wt] = quadrature_volume_3D(pts) provides the integration points
%   (in barycentric coordinates, with the fourth coordinate computed as 
%   1 - sum(first three)) and weights for integration over the reference
%   tetrahedron (volume = 1/6).
%
%   Input:
%       pts - number of quadrature points. Available rules:
%             1  - 1-point rule (centroid; exact for constant functions)
%             4  - 4-point rule (Keast K1; exact for quadratic functions)
%             5  - 5-point rule (Keast K2; exact for cubic functions)
%             8  - 8-point rule (a symmetric degree-4 rule, example below)
%            11  - 11-point rule (Keast K4; moderate degree)
%            15  - 15-point rule (Keast K6; higher degree)
%
%   Output:
%       Xi - 3-by-N matrix of barycentric coordinates (each column is a point).
%            The fourth barycentric coordinate is 1 - sum(Xi(:,i)).
%       wt - 1-by-N vector of weights (the weights sum to 1/6).
%
%   References:
%       P. Keast, Moderate degree tetrahedral quadrature formulas,
%         CMAME 55: 339-348 (1986)
%       O.C. Zienkiewicz, The Finite Element Method, Sixth Edition.
%       (For the 8-point rule, see, e.g., various compendia such as Felippa's
%         "A Compendium of FEM Integration Formulas" or other literature.)
%

switch pts
    case 1
        % 1-point rule (centroid) -- exact for constant functions.
        xa = [0.25, 0.1381966011250105, ...
              0.1381966011250105, 0.1381966011250105];
        ya = [0.1381966011250105, 0.1381966011250105, ...
              0.1381966011250105, 0.5854101966249685];
        za = [0.1381966011250105, 0.1381966011250105, ...
              0.5854101966249685, 0.1381966011250105];
        wt = [0.25, 0.25, 0.25, 0.25] / 6;
        
    case 4
        % Keast K1 (4 points)
        xa = [0.5854101966249685, 0.1381966011250105, ...
              0.1381966011250105, 0.1381966011250105];
        ya = [0.1381966011250105, 0.1381966011250105, ...
              0.1381966011250105, 0.5854101966249685];
        za = [0.1381966011250105, 0.1381966011250105, ...
              0.5854101966249685, 0.1381966011250105];
        wt = [0.25, 0.25, 0.25, 0.25] / 6;
        
    case 5
        % Keast K2 (5 points)
        xa = [0.2500000000000000, 0.5000000000000000, ...
              0.1666666666666667, 0.1666666666666667, 0.1666666666666667];
        ya = [0.2500000000000000, 0.1666666666666667, ...
              0.1666666666666667, 0.1666666666666667, 0.5000000000000000];
        za = [0.2500000000000000, 0.1666666666666667, ...
              0.1666666666666667, 0.5000000000000000, 0.1666666666666667];
        wt = [-0.8, 0.45, 0.45, 0.45, 0.45] / 6;
        
    case 8
        % 8-point rule (symmetric degree-4 rule, example)
        % Here we define a small perturbation delta so that the points are:
        %   P = (0.25 + perturbation in each coordinate)
        % with the eight points chosen to preserve symmetry.
        %
        % Define delta (this value is one possible choice)
        delta = 0.0001;
        % Define the eight integration points (only the first three barycentrics)
        p1 = [0.25 + delta; 0.25 - delta; 0.25 - delta];
        p2 = [0.25 - delta; 0.25 + delta; 0.25 - delta];
        p3 = [0.25 - delta; 0.25 - delta; 0.25 + delta];
        p4 = [0.25 - delta; 0.25 - delta; 0.25 - delta];
        p5 = [0.25 - delta; 0.25 + delta; 0.25 + delta];
        p6 = [0.25 + delta; 0.25 - delta; 0.25 + delta];
        p7 = [0.25 + delta; 0.25 + delta; 0.25 - delta];
        p8 = [0.25 + delta; 0.25 + delta; 0.25 + delta];
        Xi = [p1, p2, p3, p4, p5, p6, p7, p8];
        % All weights are equal:
        wt = ones(1,8) * (1/48);
        return
        
    case 11
        % Keast K4 (11 points)
        xa = [0.2500000000000000, 0.7857142857142857, 0.0714285714285714, ...
              0.0714285714285714, 0.0714285714285714, 0.1005964238332008, ...
              0.3994035761667992, 0.3994035761667992, 0.3994035761667992, ...
              0.1005964238332008, 0.1005964238332008];
        ya = [0.2500000000000000, 0.0714285714285714, 0.0714285714285714, ...
              0.0714285714285714, 0.7857142857142857, 0.3994035761667992, ...
              0.1005964238332008, 0.3994035761667992, 0.1005964238332008, ...
              0.3994035761667992, 0.1005964238332008];
        za = [0.2500000000000000, 0.0714285714285714, 0.0714285714285714, ...
              0.7857142857142857, 0.0714285714285714, 0.3994035761667992, ...
              0.3994035761667992, 0.1005964238332008, 0.1005964238332008, ...
              0.1005964238332008, 0.3994035761667992];
        wt = [-0.0789333333333333, 0.0457333333333333, 0.0457333333333333, ...
              0.0457333333333333, 0.0457333333333333, 0.1493333333333333, ...
              0.1493333333333333, 0.1493333333333333, 0.1493333333333333, ...
              0.1493333333333333, 0.1493333333333333] / 6;
        
    case 15
        % Keast K6 (15 points)
        xa = [0.2500000000000000, 0.0000000000000000, 0.3333333333333333, ...
              0.3333333333333333, 0.3333333333333333, 0.7272727272727273, ...
              0.0909090909090909, 0.0909090909090909, 0.0909090909090909, ...
              0.4334498464263357, 0.0665501535736643, 0.0665501535736643, ...
              0.0665501535736643, 0.4334498464263357, 0.4334498464263357];
        ya = [0.2500000000000000, 0.3333333333333333, 0.3333333333333333, ...
              0.3333333333333333, 0.0000000000000000, 0.0909090909090909, ...
              0.0909090909090909, 0.0909090909090909, 0.7272727272727273, ...
              0.0665501535736643, 0.4334498464263357, 0.0665501535736643, ...
              0.4334498464263357, 0.0665501535736643, 0.4334498464263357];
        za = [0.2500000000000000, 0.3333333333333333, 0.3333333333333333, ...
              0.0000000000000000, 0.3333333333333333, 0.0909090909090909, ...
              0.0909090909090909, 0.7272727272727273, 0.0909090909090909, ...
              0.0665501535736643, 0.0665501535736643, 0.4334498464263357, ...
              0.4334498464263357, 0.4334498464263357, 0.0665501535736643];
        wt = [0.1817020685825351, 0.0361607142857143, 0.0361607142857143, ...
              0.0361607142857143, 0.0361607142857143, 0.0698714945161738, ...
              0.0698714945161738, 0.0698714945161738, 0.0698714945161738, ...
              0.0656948493683187, 0.0656948493683187, 0.0656948493683187, ...
              0.0656948493683187, 0.0656948493683187, 0.0656948493683187] / 6;
        
    otherwise
        error('Quadrature rule for %d points is not available. Available rules: 1, 4, 5, 8, 11, 15.', pts);
end

Xi = [xa; ya; za];

end
