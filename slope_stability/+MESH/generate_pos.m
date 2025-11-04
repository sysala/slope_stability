function generate_pos(elem,coord,Xi,B,U,WEIGHT,file_name,L1_elem_size,size_mult,max_h)

%U3, coord, elem, B

iota = [1; 1; 1; 0; 0; 0];
VOL = iota * iota';
DEV = diag([1, 1, 1, 1/2, 1/2, 1/2]) - VOL / 3;

% Compute the strain tensor at all integration points.
E = B * U(:);
E = reshape(E, 6, []);  % Each column corresponds to an integration point.

% Compute the deviatoric part and its norm.
dev_E = DEV * E;
norm_E = sqrt(max(0, sum(E .* dev_E)));

[transformed_points] = ASSEMBLY.integration_points(elem,coord,Xi);


unscaled_diam = (1 ./ (norm_E)).^(1/3);
domain_volume = sum(WEIGHT);


ratio = (sum(WEIGHT./unscaled_diam.^3)/(domain_volume/(L1_elem_size/size_mult)))^(1/3);

unscaled_diam = unscaled_diam * ratio;
unscaled_diam(unscaled_diam>max_h)=max_h;
ratio = (sum(WEIGHT./unscaled_diam.^3)/(domain_volume/(L1_elem_size/size_mult)))^(1/3);
unscaled_diam = unscaled_diam * ratio;

MESH.create_pos_file(transformed_points, unscaled_diam, ['pos_files/' file_name '.pos'])

end