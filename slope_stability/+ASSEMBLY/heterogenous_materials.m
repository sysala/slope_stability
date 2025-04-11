function [c0, phi, psi, shear, bulk, lame, gamma] = heterogenous_materials(mat_identifier, n_q, materials)
%--------------------------------------------------------------------------
% heterogenous_materials assigns material properties for heterogeneous media.
%
% This function initializes and assigns material properties for each element
% based on a material identifier array and a cell array of material property
% structures. The properties include c0, phi, psi, shear modulus, bulk modulus,
% Lame's coefficient (lambda), and gamma. Angular values (phi and psi) are
% converted from degrees to radians. The properties are then expanded for
% quadrature by repeating each value n_q times.
%
% INPUT ARGUMENTS:
%   mat_identifier - Array of material identifiers for each element. It is
%                    assumed that these identifiers start from 0.
%   n_q            - Number of quadrature points per element.
%   materials      - Cell array of structures, where each structure contains 
%                    material properties for a specific material. Each structure
%                    must include fields: c0, phi, psi, young, poisson, gamma.
%
% OUTPUT:
%   c0     - Array of c0 values for each integration point, size: (1, n_int*n_q)
%   phi    - Array of phi values (in radians) for each integration point.
%   psi    - Array of psi values (in radians) for each integration point.
%   shear  - Array of shear modulus values for each integration point.
%   bulk   - Array of bulk modulus values for each integration point.
%   lame   - Array of Lame's coefficient (lambda) values for each integration point.
%   gamma  - Array of gamma values for each integration point.
%
%--------------------------------------------------------------------------

% Number of elements (integration points before quadrature expansion)
n_int = length(mat_identifier);

% Initialize material property arrays (one entry per element)
c0    = zeros(1, n_int);
phi   = zeros(1, n_int);
psi   = zeros(1, n_int);
shear = zeros(1, n_int);
bulk  = zeros(1, n_int);
lame  = zeros(1, n_int);
gamma = zeros(1, n_int);

% Assign properties to each element based on the material identifier
for i = 1:n_int
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
    gamma(i) = materials{mat_id}.gamma;
end

% Expand the material property arrays to match the number of quadrature points.
c0    = kron(c0, ones(1, n_q));
phi   = kron(phi, ones(1, n_q));
psi   = kron(psi, ones(1, n_q));
shear = kron(shear, ones(1, n_q));
bulk  = kron(bulk, ones(1, n_q));
lame  = kron(lame, ones(1, n_q));
gamma = kron(gamma, ones(1, n_q));
end
