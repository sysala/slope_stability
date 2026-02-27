function [Q_w, pw_D] = seepage_boundary_3D_hetero_comsol(coord, surf, triangle_labels, grho)

%--------------------------------------------------------------------------
% This function specifies homogeneous and nonhomogeneous Dirichlet boundary
% conditions for the unconfined seepage problem corresponding to the
% heterogeneous concave slope. This function is not universal. It is
% prepared only for a particular problem!
%
% INPUT ARGUMENTS:
%   coord: Coordinates of degrees of freedom, computed from mesh nodes and element_order.
%   surf: Surface connectivity matrix (triangle2node), each row contains dof indices for a surface.
%         surface in question is mostly boundary but also not boundary
%   triangle_labels: An array specifying a surface type for each surface
%                    triangle. The coordinate y is the height of the body.
%   grho: specific weight of water in kPa, i.e., grho=9.81
%
% OUTPUT:
%   Q_darcy: Logical array excluding nodes with Dirichlet conditions for
%            pore pressure
%   pw_D: Nonhomogeneous part of the pressure (problem dependent)
%
%--------------------------------------------------------------------------

% number of nodes
n_n = size(coord,2);

% water levels prescribed on opposite sides of the slope
y_free_water_level=35;   % lower water level (on the slope part)
y_porous_water_level=55; % higher water level (on the part opposite to the slope)

% nonhomogeneous part of the pressure - initialization
pw_D=zeros(1,n_n);
x = coord(1, :);
y = coord(2, :);
z = coord(3, :);
boundary_nodes = unique(surf(:));

% zero Dirichlet bc on top face, no changes in pw_D needed
Q_dry=false(1,n_n);
tmp = surf(:, triangle_labels == 6); % face with y=y_max - top of the slope
Q_dry(tmp(:)) = 1;

% Non-homogeneous Dirichlet boundary conditions on the wet part of the face with x == x_max
Q_wet1=false(1,n_n);
tmp = surf(:, triangle_labels == 2); % x == x_max
nodes_tmp = unique(tmp(:));
selected_nodes_wet = nodes_tmp( y(nodes_tmp) < y_porous_water_level );
Q_wet1(selected_nodes_wet) = 1;
pw_D(Q_wet1)=grho*(y_porous_water_level - coord(2,Q_wet1));
% Homogeneous Dirichlet boundary condition on the dry part of x == x_max
dry = nodes_tmp( y(nodes_tmp) >= y_porous_water_level );
Q_dry(dry) = 1;

% Non-homogeneous Dirichlet boundary conditions on the remaining faces
Q_wet2=false(1,n_n);
% find the triangles living on the convex slopes
% triangles as vertices
triangles = surf(1:3, :);
% vertex coordinates
v1 = coord(:, triangles(1, :));
v2 = coord(:, triangles(2, :));
v3 = coord(:, triangles(3, :));
% edges as "vectors"
e1 = v2 - v1;
e2 = v3 - v1;
normals = cross(e1, e2, 1);
% wet part of the slope face
tol = 1e-1; % problem dependant, brittle
condition = all(abs(normals) > tol, 1);
selected_triangles = surf(:, condition);
nodes_tmp = unique(selected_triangles(:));
% now, nodes_tmp contain all nodes that live in the two parallel planes
% one is the outer boundary, another is the inner boundary of the cover
tol = 1e-6;
C = [55, 30, 0];
T = [115, 60, 0];
A_left = [30, 30, 43.3];
A_right = A_left;
A_right(3) = -A_right(3);
normal_left = cross(T - C, A_left - C);
normal_right = cross(T - C, A_right - C);
X = coord(:, nodes_tmp).';
V = X - C;
d_left = abs(V * normal_left.');
d_right = abs(V * normal_right.');
mask = (d_left < tol) | (d_right < tol);

nodes_tmp = nodes_tmp(mask);

% assign the actual condition
selected_nodes = nodes_tmp( y(nodes_tmp) < y_free_water_level );
Q_wet2(selected_nodes) = 1;
% dry part of the slope face
selected_nodes_dry = nodes_tmp( y(nodes_tmp) >= y_free_water_level );
Q_dry(selected_nodes_dry) = 1;

% water bed
tol = 1e-1;  % problem dependant, brittle
y_bed = 30;
% needs fix, also includes the inner boundary
selected = boundary_nodes( abs(y(boundary_nodes) - y_bed) < tol);
% selected now contains all nodes with y = y_bed, now to reduce them to
% outer nodes only
% use geometry info because fuck abstraction and separation
tol = 1e-10;
C = [55, 30, 0];
A_left = [30, 30, 43.3];
A_right = A_left;
A_right(3) = -A_right(3);
n = [0, 1, 0];

X = coord(:, selected).';

dL = A_left  - C;     % left boundary direction
dR = A_right - C;     % right boundary direction

VL = X - C;           % Nx3

kL = cross(n, dL);    % 1x3
kR = cross(n, dR);    % 1x3
kL = kL/norm(kL);
kR = kR/norm(kR);

sL = VL * kL.';       % Nx1
sR = VL * kR.';
sL = sL/norm(sL);
sR = sR/norm(sR);

mask = (sL < tol) & (sR > tol);

selected = selected(mask);

Q_wet2(selected) = 1;

% the front x == 0 face
tmp = surf(:, triangle_labels == 1);  % face with x=x_max
Q_wet2(tmp(:)) = 1;
pw_D(Q_wet2)=grho*(y_free_water_level - coord(2,Q_wet2));


% Bottom face -> zero Neumann
% Side faces -> zero Neumann

% The logical array Q_w
Q_w=true(1,n_n);
Q_w(Q_dry) = 0;
Q_w(Q_wet2) = 0;
Q_w(Q_wet1) = 0;
