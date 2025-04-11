function f_V = vector_volume_3D(ELEM, COORD, f_V_int, HatP, WEIGHT)
%--------------------------------------------------------------------------
% vector_volume_3D assembles the vector of volume forces for a 3D domain.
%
% This function computes the nodal volume force vector f_V by integrating 
% the force contributions from all elements. The force values at the integration
% points are distributed to the nodes using the basis functions evaluated 
% at the quadrature points.
%
% The volume force vector f_V is assembled in a sparse format with size (3, n_n),
% where n_n is the number of nodes.
%
% INPUT ARGUMENTS:
%   ELEM    - Element connectivity matrix, size: (n_p, n_e)
%             n_e: number of elements, n_p: number of nodes per element.
%   COORD   - Nodal coordinates, size: (3, n_n)
%   f_V_int - Volume force values at integration points, size: (3, n_int)
%             n_int = n_e * n_q, where n_q is the number of quadrature points.
%   HatP    - Basis functions evaluated at quadrature points, size: (n_p, n_q)
%   WEIGHT  - Integration weight coefficients, size: (1, n_int)
%
% OUTPUT:
%   f_V     - Assembled vector of volume forces, size: (3, n_n)
%
%--------------------------------------------------------------------------

    %-------------------------------------------------------------------------
    % Auxiliary notation
    %-------------------------------------------------------------------------
    n_n = size(COORD, 2);   % Number of nodes (including midpoints)
    n_e = size(ELEM, 2);    % Number of elements
    n_p = size(ELEM, 1);    % Number of vertices per element
    n_q = size(HatP, 2);    % Number of quadrature points per element
    n_int = n_e * n_q;      % Total number of integration points

    % Extension of the input array HatP by replication to match the total 
    % number of integration points.
    % Resulting size of HatPhi is (n_p, n_int)
    HatPhi = repmat(HatP, 1, n_e);

    %-------------------------------------------------------------------------
    % Assemble the vector of volume forces (f_V) with size (3, n_n)
    %-------------------------------------------------------------------------
    % Compute contributions at integration points for each force component.
    % Each vF* has size (n_p, n_int)
    vF1 = HatPhi .* (ones(n_p, 1) * (WEIGHT .* f_V_int(1, :)));
    vF2 = HatPhi .* (ones(n_p, 1) * (WEIGHT .* f_V_int(2, :)));
    vF3 = HatPhi .* (ones(n_p, 1) * (WEIGHT .* f_V_int(3, :)));

    % Determine the row indices (iF) for sparse assembly.
    % iF is a matrix of ones with size (n_p, n_int)
    iF = ones(n_p, n_int);

    % Determine the column indices (jF) for sparse assembly.
    % jF is obtained by replicating the connectivity matrix for each quadrature point.
    jF = kron(ELEM, ones(1, n_q));

    % Assemble the volume force vector using the sparse command.
    % Duplicate entries are automatically summed.
    f_V = [ sparse(iF(:), jF(:), vF1(:), 1, n_n);
            sparse(iF(:), jF(:), vF2(:), 1, n_n);
            sparse(iF(:), jF(:), vF3(:), 1, n_n) ];
end
