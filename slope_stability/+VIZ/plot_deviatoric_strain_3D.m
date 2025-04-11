function fig = plot_deviatoric_strain_3D(U, coord, elem, B)
%--------------------------------------------------------------------------
% plot_deviatoric_strain_3D visualizes the deviatoric strain field on a 3D mesh.
%
% This function computes the element-wise deviatoric strain from the displacement
% field U and the strain-displacement matrix B using the function 
% VIZ.get_elem_stress_3D. It then extracts the boundary faces and associated values 
% for the computed deviatoric strain via VIZ.boundary_stress_from_elements_3D.
%
% Finally, the function calls VIZ.draw_quantity_3D to render a 3D plot, where the 
% color represents the magnitude of the deviatoric strain.
%
% INPUTS:
%   U     - Displacement field, size (3, n_n), where n_n is the number of nodes.
%   coord - Nodal coordinates, size (3, n_n).
%   elem  - Element connectivity matrix.
%   B     - Strain-displacement matrix used to compute the strain.
%
% OUTPUT:
%   fig   - Handle to the generated figure.
%
%--------------------------------------------------------------------------
%%
% Compute element-wise deviatoric strain values.
[elem_values] = VIZ.get_elem_stress_3D(U, B);
%%
% Extract boundary faces, boundary node coordinates, and corresponding strain values.
[boundary_faces, coord_boundary, values_boundary] = VIZ.boundary_stress_from_elements_3D(coord, elem, elem_values);
%%
% Plot the deviatoric strain field.
% Here, the third input (quantity to be plotted on the faces) is set to 0,
% so that only the scalar values (values_boundary) determine the color.
fig = VIZ.draw_quantity_3D(coord_boundary, boundary_faces, 0 * values_boundary, values_boundary);
title("Deviatoric strain")
end
