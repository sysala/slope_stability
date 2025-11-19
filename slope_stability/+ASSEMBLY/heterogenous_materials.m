function [c0, phi, psi, shear, bulk, lame, gamma] = ...
         heterogenous_materials(mat_identifier, saturation, n_q, materials)
%--------------------------------------------------------------------------
% heterogenous_materials assigns material properties for heterogeneous media.
%
% This function initializes and assigns material properties for each element
% based on a material identifier array and a cell array of material property
% structures. The properties include c0, phi, psi, shear modulus, bulk modulus,
% Lame's coefficient (lambda), gamma_sat and gamma_unsat. Angular values 
% (phi and psi) are converted from degrees to radians. The properties are 
% then expanded for quadrature by repeating each value n_q times. The specific
% weight gamma at a particular integration point is determined from 
% gamma_sat, gamma_unsat, and a given logical array saturation.
%
% INPUT ARGUMENTS:
%   mat_identifier - Array of material identifiers for each element. It is
%                    assumed that these identifiers start from 0.
%   saturation     - a logical array indicating integration points 
%                    where the body is saturated.
%   n_q            - Number of quadrature points per element.
%   materials      - Cell array of structures, where each structure contains 
%                    material properties for a specific material. Each structure
%                    must include fields: c0, phi, psi, young, poisson, gamma.
%
% OUTPUT:
%   c0     - Array of c0 values for each integration point, size: (1, n_e*n_q)
%   phi    - Array of phi values (in radians) for each integration point.
%   psi    - Array of psi values (in radians) for each integration point.
%   shear  - Array of shear modulus values for each integration point.
%   bulk   - Array of bulk modulus values for each integration point.
%   lame   - Array of Lame's coefficient (lambda) values for each integration point.
%   gamma  - Array of specific weights for each integration point.
%
%--------------------------------------------------------------------------

% Numbers of elements and integration points
n_e = length(mat_identifier);
n_int=n_e*n_q;

% Initialize material property arrays (one entry per element)
c0    = zeros(1, n_e);
phi   = zeros(1, n_e);
psi   = zeros(1, n_e);
shear = zeros(1, n_e);
bulk  = zeros(1, n_e);
lame  = zeros(1, n_e);
gamma_sat = zeros(1, n_e);
gamma_unsat = zeros(1, n_e);

% Assign properties to each element based on the material identifier
for i = 1:n_e
    % Adjust index since mat_identifier is assumed to start from 0
    mat_id = mat_identifier(i) + 1;
    
    % Retrieve material properties from the corresponding structure.
    c0(i)   = materials{mat_id}.c0;
    % Convert angles from degrees to radians.
    phi(i)  = materials{mat_id}.phi * pi / 180;
    psi(i)  = materials{mat_id}.psi * pi / 180;
    
    % Calculate elastic properties using Young's modulus and Poisson's ratio.
    young   = materials{mat_id}.young;
    poisson = materials{mat_id}.poisson;
    shear(i) = young / (2 * (1 + poisson));      % Shear modulus
    bulk(i)  = young / (3 * (1 - 2 * poisson));    % Bulk modulus
    lame(i)  = bulk(i) - 2 * shear(i) / 3;          % Lame's first parameter (lambda)
    
    % Assign gamma directly.
    gamma_sat(i) = materials{mat_id}.gamma_sat;
    gamma_unsat(i) = materials{mat_id}.gamma_unsat;
end

% Expand the material property arrays to match the number of quadrature points.
c0    = kron(c0, ones(1, n_q));
phi   = kron(phi, ones(1, n_q));
psi   = kron(psi, ones(1, n_q));
shear = kron(shear, ones(1, n_q));
bulk  = kron(bulk, ones(1, n_q));
lame  = kron(lame, ones(1, n_q));
gamma_sat = kron(gamma_sat, ones(1, n_q));
gamma_unsat = kron(gamma_unsat, ones(1, n_q));

% Determination of the array with the specific weight
gamma = zeros(1, n_int);
gamma(saturation)=gamma_sat(saturation);
gamma(~saturation)=gamma_unsat(~saturation);
end
