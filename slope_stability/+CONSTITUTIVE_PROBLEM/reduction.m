function [c_bar, sin_phi] = reduction(c0, phi, psi, lambda, Davis_type)
%--------------------------------------------------------------------------
% reduction reduces strength parameters according to a strength reduction factor
% and the prescribed Davis' approach.
%
% This function computes reduced strength parameters based on the input
% effective cohesion (c0), effective friction angle (phi), dilatancy angle (psi),
% and the strength reduction factor (lambda) following one of three Davis'
% approaches ('A', 'B', or 'C'). The reduction results in a modified cohesion
% and friction angle, which are then used to compute c_bar and sin_phi.
%
% INPUT ARGUMENTS:
%   c0         - Effective cohesion at integration points.
%   phi        - Effective friction angle at integration points (in radians).
%   psi        - Dilatancy angle at integration points (in radians).
%   lambda     - Strength reduction factor.
%   Davis_type - Character indicating the Davis' approach to be used.
%                Supported options: 'A', 'B', 'C'.
%
% OUTPUT:
%   c_bar  - Reduced cohesion parameter, computed as 2*c0_lambda*cos(phi_lambda).
%   sin_phi- Sine of the reduced friction angle, sin(phi_lambda).
%
%--------------------------------------------------------------------------
  switch(Davis_type)
    case 'A'
        % Approach A: Direct reduction based on original angles.
        beta = cos(phi) .* cos(psi) ./ (1 - sin(phi) .* sin(psi));
        c0_lambda = beta .* c0 / lambda;
        phi_lambda = atan(beta .* tan(phi) / lambda);
        c_bar = 2 * c0_lambda .* cos(phi_lambda);   
        sin_phi = sin(phi_lambda); 
    case 'B'
        % Approach B: Reduction applied to the angles first.
        c01 = c0 / lambda;
        phi1 = atan(tan(phi) / lambda);
        psi1 = atan(tan(psi) / lambda);
        beta = cos(phi1) .* cos(psi1) ./ (1 - sin(phi1) .* sin(psi1));
        c0_lambda = beta .* c01;
        phi_lambda = atan(beta .* tan(phi1));
        c_bar = 2 * c0_lambda .* cos(phi_lambda);   
        sin_phi = sin(phi_lambda); 
    case 'C'
        % Approach C: Conditional beta based on reduced friction angle versus dilatancy.
        c01 = c0 / lambda;
        phi1 = atan(tan(phi) / lambda);
        if phi1 > psi
            beta = cos(phi1) .* cos(psi) ./ (1 - sin(phi1) .* sin(psi));
        else
            beta = 1;
        end
        c0_lambda = beta .* c01;
        phi_lambda = atan(beta .* tan(phi1));
        c_bar = 2 * c0_lambda .* cos(phi_lambda);   
        sin_phi = sin(phi_lambda); 
    otherwise
        disp('Incorrect choice of the Davis approach.');
  end        
  
end
