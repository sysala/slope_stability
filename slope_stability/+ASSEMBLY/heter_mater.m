function [c0, phi, psi, gamma_sat, gamma_unsat] = heter_mater(mater, n_q, materials)
%--------------------------------------------------------------------------
% heter_mater assigns material properties for heterogeneous media.
%
% This function initializes and assigns material properties for each element
% based on a material identifier array and a cell array of material property
% structures. The properties include c0, phi, psi, gamma_sat, and gamma_unsat.
% The properties are then expanded for quadrature by repeating each value n_q times.
%
% INPUT ARGUMENTS:
%   mater      - Array of material identifiers for each element.
%   n_q        - Number of quadrature points per element.
%   materials  - Cell array of structures, where each structure contains
%                material properties for a specific material.
%
% OUTPUT:
%   c0         - Array of c0 values for each integration point.
%   phi        - Array of phi values (in radians) for each integration point.
%   psi        - Array of psi values (in radians) for each integration point.
%   gamma_sat  - Array of saturated unit weights for each integration point.
%   gamma_unsat- Array of unsaturated unit weights for each integration point.
%
%--------------------------------------------------------------------------

% Number of elements (integration points before quadrature expansion)
n_e = length(mater);

% Initialize material property arrays (one entry per element)
c0         = zeros(1, n_e);
phi        = zeros(1, n_e);
psi        = zeros(1, n_e);
gamma_sat  = zeros(1, n_e);
gamma_unsat = zeros(1, n_e);

% Assign properties to each element based on the material identifier
for i = 1:n_e
    mat_id = mater(i); % Material identifier
    mat = materials{mat_id};
    
    c0(i)         = mat.c0;
    phi(i)        = mat.phi; % Already in radians
    psi(i)        = mat.psi; % Already in radians
    gamma_sat(i)  = mat.gamma_sat;
    gamma_unsat(i) = mat.gamma_unsat;
end

% Expand the material property arrays to match the number of quadrature points.
c0         = kron(c0, ones(1, n_q));
phi        = kron(phi, ones(1, n_q));
psi        = kron(psi, ones(1, n_q));
gamma_sat  = kron(gamma_sat, ones(1, n_q));
gamma_unsat = kron(gamma_unsat, ones(1, n_q));

end
