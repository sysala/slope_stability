function [Q_darcy, pw_D] = darcy_boundary_3D_hetero(coord, surf, triangle_labels, grho)
% Load mesh data from a GMSH file with water levels.

% OUTPUTS:
%   coord: Coordinates of degrees of freedom, computed from mesh nodes and element_order.
%   elem: Element connectivity matrix (tetra2node), each row contains dof indices for an element.
%   surf: Surface connectivity matrix (triangle2node), each row contains dof indices for a surface.
%   Q: Logical arrat indicating dofs restricted by boundary conditions.
%   material: Material indices for each element

n_n = size(coord,2);

z_water_height=35;
z_solid_water_level=50;


% name2tag = {'bottom_layer': 1, 'middle_layer': 2, 'top_layer': 3, 'cover': 4,
%             'dry_slope': 5, 'wet_slope': 6,
%             'dry_solid': 7, 'wet_solid': 8,
%             'x_max': 9, 'y0': 10, 'ymax': 11, 'z0': 12, 'zmax': 13, 'zwater_bed': 14}

pw_D=zeros(1,n_n);

% Dirichlet boundary conditions for pressure (problem dependent)
Q_open=false(1,n_n);

tmp = surf(:, triangle_labels == 13); %zmax
Q_open(tmp(:)) = 1;
tmp = surf(:, triangle_labels == 5); %dry slope
Q_open(tmp(:)) = 1;
tmp = surf(:, triangle_labels == 7); %dry solid
Q_open(tmp(:)) = 1;

Q_solid_wet=false(1,n_n);
tmp = surf(:, triangle_labels == 8); %wet solid
Q_solid_wet(tmp(:)) = 1;
pw_D(Q_solid_wet)=grho*(z_solid_water_level - coord(2,Q_solid_wet));

Q_open_wet=false(1,n_n);
tmp = surf(:, triangle_labels == 6); %wet slope
Q_open_wet(tmp(:)) = 1;
tmp = surf(:, triangle_labels == 14); %water bed
Q_open_wet(tmp(:)) = 1;
tmp = surf(:, triangle_labels == 9); %x max
Q_open_wet(tmp(:)) = 1;
pw_D(Q_open_wet)=grho*(z_water_height - coord(2,Q_open_wet));

Q_darcy=true(1,n_n);
Q_darcy(Q_open) = 0;
Q_darcy(Q_open_wet) = 0;
Q_darcy(Q_solid_wet) = 0;