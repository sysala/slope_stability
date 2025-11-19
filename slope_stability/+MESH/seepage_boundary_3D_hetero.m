function [Q_w, pw_D] = seepage_boundary_3D_hetero(coord, surf, triangle_labels, grho)

%--------------------------------------------------------------------------
% This function specifies homogeneous and nonhomogeneous Dirichlet boundary
% conditions for the unconfined seepage problem corresponding to the
% heterogeneous concave slope. This function is not universal. It is
% prepared only for a particular problem!
%
% INPUT ARGUMENTS:
%   coord: Coordinates of degrees of freedom, computed from mesh nodes and element_order.
%   surf: Surface connectivity matrix (triangle2node), each row contains dof indices for a surface.
%   triangle_labels: An array specifying a surface type for each surface
%                    triangle. The coordinate y is the height of the body.
%      1 - bottom layer
%      2 - middle layer
%      3 - top layer
%      4 - cover
%      5 - dry part of the slope face
%      6 - wet part of the slope face
%      7 - dry part of the face with x=0
%      8 - wet part of the face with x=0
%      9 - face with x=x_max
%     10 - face with z=0
%     11 - face with z=z_max
%     12 - face with y=0 - bottom of the foundation
%     13 - face with y=y_max - top of the slope
%     14 - water bed face (constant y representing height of the foundation) 
%   grho: specific weight of water in kPa, i.e., grho=9.81
%
% OUTPUT:
%   Q_darcy: Logical array excluding nodes with Dirichlet conditions for
%            pore pressures
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

% Homogeneous Dirichlet boundary conditions for dry part of the boundary
Q_dry=false(1,n_n);
tmp = surf(:, triangle_labels == 13); % face with y=y_max - top of the slope
Q_dry(tmp(:)) = 1;
tmp = surf(:, triangle_labels == 5);  % dry part of the slope face
Q_dry(tmp(:)) = 1;
tmp = surf(:, triangle_labels == 7);  % dry part of the face with x=0
Q_dry(tmp(:)) = 1;

% Non-homogeneous Dirichlet boundary conditions on the wet part of the face with x=0
Q_wet1=false(1,n_n);
tmp = surf(:, triangle_labels == 8); % wet part of the face with x=0
Q_wet1(tmp(:)) = 1;
pw_D(Q_wet1)=grho*(y_porous_water_level - coord(2,Q_wet1));

% Non-homogeneous Dirichlet boundary conditions on the remaining wet faces
Q_wet2=false(1,n_n);
tmp = surf(:, triangle_labels == 6);  % wet part of the slope face
Q_wet2(tmp(:)) = 1;
tmp = surf(:, triangle_labels == 14); % water bed face
Q_wet2(tmp(:)) = 1;
tmp = surf(:, triangle_labels == 9);  % face with x=x_max
Q_wet2(tmp(:)) = 1;
pw_D(Q_wet2)=grho*(y_free_water_level - coord(2,Q_wet2));

% The logical array Q_w
Q_w=true(1,n_n);
Q_w(Q_dry) = 0;
Q_w(Q_wet2) = 0;
Q_w(Q_wet1) = 0;