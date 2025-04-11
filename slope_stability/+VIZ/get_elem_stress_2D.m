function [elem_values] = get_elem_stress_2D(U, B, elem)
%--------------------------------------------------------------------------
% get_elem_stress_3D computes the deviatoric strain norm at integration points
% for a 3D finite element.
%
% This function calculates the norm of the deviatoric part of the strain tensor
% at each integration point. The strain tensor E is computed as:
%
%       E = B * U(:),
%
% where U is the displacement field and B is the strain-displacement matrix.
%
% The deviatoric part of the strain is obtained by applying the projection operator:
%
%       DEV = I_dev = diag([1, 1, 1, 1/2, 1/2, 1/2]) - (1/3) * (iota * iota'),
%
% with iota = [1; 1; 1; 0; 0; 0]. The norm of the deviatoric strain is then
% computed in an A-inner product sense (here simply as sqrt(sum(E .* (DEV*E)))),
% and subsequently reshaped to form an array of element-wise values.
%
% Finally, the rows of the reshaped matrix are reordered according to a fixed
% pattern.
%
% INPUTS:
%   U           - Displacement field (vector or matrix form).
%   B           - Strain-displacement matrix.
%   elem        - Element connectivity matrix, size (3 x n), where n is the number of elements.
%
% OUTPUT:
%   elem_values - Matrix containing the deviatoric strain norms for each element.
%                 The output is reshaped and reordered based on the integration scheme.
%
%--------------------------------------------------------------------------
iota = [1;1;0;1];
VOL = iota * iota';
DEV = diag([1,1,1/2,1])-VOL/3;

% Compute the strain tensor at all integration points.
E = B * U(:);   % strain tensors at integration points
E = reshape(E, 3, []);  % Each column corresponds to an integration point.
E(end+1,:)=0;   % trial strain

% Compute the deviatoric part and its norm.
dev_E = DEV * E;
norm_E = sqrt(max(0, sum(E .* dev_E)));

% Reshape the norm vector into a matrix with 11 rows (integration points per element).
if size(elem,1) == 6
    elem_values = reshape(norm_E, 7, []);
elseif size(elem,1) == 15
    elem_values = reshape(norm_E, 12, []);
else
    elem_values = repmat(reshape(norm_E, 1, []),3,1);
end
% Reorder the rows to a specific pattern.
if size(elem,1) == 6
    elem_values = elem_values([1 2 3 5 6 4], :);
elseif size(elem,1) == 15
    elem_values_ = elem_values([1 2 3], :);
    elem_values_(4,:) = sum(elem_values([9,10],:))/2;
    elem_values_(5,:) = sum(elem_values([8,11],:))/2;
    elem_values_(6,:) = sum(elem_values([7,12],:))/2;
    elem_values = elem_values_;
end

end
