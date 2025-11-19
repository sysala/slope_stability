function [coord, elem, surf, Q, material, boundary] = load_mesh_P2(file_path, boundary_type)
%--------------------------------------------------------------------------
% load_mesh_P2 loads mesh data for quadratic (P2) finite elements.
%
% This function reads mesh datasets from an HDF5 file specified by file_path.
% It loads the nodal coordinates, element connectivity, face connectivity,
% boundary conditions, and material indices. The element and face connectivity
% arrays are shifted by 1 to account for MATLAB's 1-indexing.
%
% OUTPUTS:
%   coord    - Nodal coordinates (3 x n_nodes), reordered as [x; z; y].
%   elem     - Element connectivity matrix for P2 elements, as double.
%   surf     - Face connectivity matrix, as double.
%   Q        - Logical array (3 x n_nodes) indicating restricted DOFs 
%              due to boundary conditions.
%   material - Material indices for each element.
%--------------------------------------------------------------------------
%

if nargin < 2
    boundary_type = 0;
end

% Load datasets from the HDF5 file.
boundary = h5read(file_path, '/boundary');
elem = h5read(file_path, '/elem') + 1;  % Adjust for MATLAB 1-indexing.
face = h5read(file_path, '/face') + 1;  % Adjust for MATLAB 1-indexing.
material = h5read(file_path, '/material');
node = h5read(file_path, '/node');

% Initialize logical array Q (3 x n_nodes) with true values.
Q = true(size(node));

% No change to first coordinate.
node(1,:) = node(1,:);

% For each boundary condition, update Q to mark the corresponding DOFs as false.
tmp = face(:, boundary == 1);
Q(1, tmp(:)) = 0;
tmp = face(:, boundary == 2);
Q(1, tmp(:)) = 0;
tmp = face(:, boundary == 3);
Q(2, tmp(:)) = 0;
tmp = face(:, boundary == 4);
Q(2, tmp(:)) = 0;
tmp = face(:, boundary == 5);
if boundary_type
    Q(:, tmp(:)) = 0;
else
    Q(3, tmp(:)) = 0;
end

% Reorder coordinates: switch order to [x; z; y] (from original [x; y; z]).
coord = double(node([1 3 2], :));
Q = Q([1 3 2], :);
elem = double(elem);
surf = double(face);
