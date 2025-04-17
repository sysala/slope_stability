function Psi = potential_3D(E, c_bar, sin_phi, shear, bulk, lame)
%--------------------------------------------------------------------------
% potential_3D computes the potential for the constitutive Drucker-Prager 
% problem using the Mohr-Coulomb yield criterion and an associative flow rule.
%
% This function constructs a potential function Psi based on the current strain
% tensor E at each integration point. The potential accounts for both elastic and 
% inelastic (plastic) responses by evaluating the invariants of the trial strain,
% determining the ordered eigenvalues, and then selecting the appropriate response 
% regime based on critical criteria.
%
% INPUT ARGUMENTS:
%   E       - Current strain tensor at integration points, size: (6, n_int)
%             The components of E are ordered as [E11; E22; E33; E12; E13; E23].
%   c_bar   - Reduced cohesion parameter at integration points, size: (1, n_int)
%   sin_phi - Sine of the reduced friction angle at integration points, size: (1, n_int)
%   shear   - Shear modulus at integration points, size: (1, n_int)
%   bulk    - Bulk modulus at integration points, size: (1, n_int)
%   lame    - Lame's coefficient (lambda) at integration points, size: (1, n_int)
%
% OUTPUT:
%   Psi     - Potential function evaluated at each integration point, 
%             size: (1, n_int)
%
%--------------------------------------------------------------------------
%
% Auxiliary definitions and trial strain evaluation
%
n_int = size(E, 2);          % Number of integration points
E_trial = E;                 % Trial strain tensor remains the same
IDENT = diag([1, 1, 1, 1/2, 1/2, 1/2]);  % Identity operator for strain conversion
E_tr = IDENT * E_trial;      % Transformed strain tensor (stress representation)

%
% Compute invariants of the trial strain tensor at each integration point.
%
I1 = E_tr(1,:) + E_tr(2,:) + E_tr(3,:);  % First invariant (trace)
I2 = E_tr(1,:) .* E_tr(2,:) + E_tr(1,:) .* E_tr(3,:) + E_tr(2,:) .* E_tr(3,:) - ...
     E_tr(4,:).^2 - E_tr(5,:).^2 - E_tr(6,:).^2;  % Second invariant
I3 = E_tr(1,:) .* E_tr(2,:) .* E_tr(3,:) - E_tr(3,:) .* E_tr(4,:).^2 - ...
     E_tr(2,:) .* E_tr(6,:).^2 - E_tr(1,:) .* E_tr(5,:).^2 + ...
     2 * E_tr(4,:) .* E_tr(5,:) .* E_tr(6,:);  % Third invariant

% Compute auxiliary quantities for Lode's angle calculation.
Q = max(0, (1/9) * (I1.^2 - 3 * I2));
R = (1/54) * (-2 * (I1.^3) + 9 * (I1 .* I2) - 27 * I3);
theta0 = zeros(1, n_int);
theta0(Q > 0) = R(Q > 0) ./ sqrt(Q(Q > 0).^3);
theta = acos(min(max(theta0, -1), 1)) / 3;  % Lode's angle

%
% Compute ordered eigenvalues of the trial strain tensors.
%
eig_1 = -2 * sqrt(Q) .* cos(theta + 2*pi/3) + I1 / 3;
eig_2 = -2 * sqrt(Q) .* cos(theta - 2*pi/3) + I1 / 3;
eig_3 = -2 * sqrt(Q) .* cos(theta) + I1 / 3;

%
% Define critical values for yield criteria.
%
% f_tr defines the yield function (elasto-plastic interface)
f_tr = 2 * shear .* ((1 + sin_phi) .* eig_1 - (1 - sin_phi) .* eig_3) + ...
       2 * (lame .* sin_phi) .* I1 - c_bar;
% Interfaces for different return types on the yield surface.
gamma_sl = (eig_1 - eig_2) ./ (1 + sin_phi);         % Smooth-left interface
gamma_sr = (eig_2 - eig_3) ./ (1 - sin_phi);         % Smooth-right interface
gamma_la = (eig_1 + eig_2 - 2 * eig_3) ./ (3 - sin_phi);  % Left-apex interface
gamma_ra = (2 * eig_1 - eig_2 - eig_3) ./ (3 + sin_phi);  % Right-apex interface

%
% Compute candidate plastic multipliers for different return mappings.
%
denom_s = 4 * lame .* sin_phi.^2 + 4 * shear .* (1 + sin_phi.^2);   
denom_l = 4 * lame .* sin_phi.^2 + shear .* (1 + sin_phi).^2 + ...
          2 * shear .* (1 - sin_phi).^2;
denom_r = 4 * lame .* sin_phi.^2 + 2 * shear .* (1 + sin_phi).^2 + ...
          shear .* (1 - sin_phi).^2;
denom_a = 4 * bulk .* sin_phi.^2;  % Denominator for apex return

lambda_s = f_tr ./ denom_s;
lambda_l = (shear .* ((1 + sin_phi) .* (eig_1 + eig_2) - 2 * (1 - sin_phi) .* eig_3) + ...
            2 * lame .* sin_phi .* I1 - c_bar) ./ denom_l;
lambda_r = (shear .* (2 * (1 + sin_phi) .* eig_1 - (1 - sin_phi) .* (eig_2 + eig_3)) + ...
            2 * lame .* sin_phi .* I1 - c_bar) ./ denom_r;
lambda_a = (2 * bulk .* sin_phi .* I1 - c_bar) ./ denom_a;

%
% Determine the potential Psi based on the response regime.
%
Psi = zeros(1, n_int);       % Initialize potential array
trace_E = eig_1 + eig_2 + eig_3;  % Trace of the trial strain

% Elastic response (no yielding)
test_el = (f_tr <= 0);
Psi(test_el) = 0.5 * lame(test_el) .* trace_E(test_el).^2 + ...
               shear(test_el) .* (eig_1(test_el).^2 + eig_2(test_el).^2 + eig_3(test_el).^2);

% Return mapping: smooth portion of the yield surface.
test_s = (lambda_s <= min(gamma_sl, gamma_sr)) & (~test_el);
Psi(test_s) = 0.5 * lame(test_s) .* trace_E(test_s).^2 + ...
              shear(test_s) .* (eig_1(test_s).^2 + eig_2(test_s).^2 + eig_3(test_s).^2) - ...
              0.5 * denom_s(test_s) .* (lambda_s(test_s).^2);
                
% Return mapping: left edge of the yield surface.
test_l = (gamma_sl < gamma_sr) & (lambda_l >= gamma_sl) & (lambda_l <= gamma_la) & ...
         ~(test_el | test_s);
Psi(test_l) = 0.5 * lame(test_l) .* trace_E(test_l).^2 + ...
              shear(test_l) .* (eig_3(test_l).^2 + 0.5 * (eig_1(test_l) + eig_2(test_l)).^2) - ...
              0.5 * denom_l(test_l) .* (lambda_l(test_l).^2);
            
% Return mapping: right edge of the yield surface.
test_r = (gamma_sl > gamma_sr) & (lambda_r >= gamma_sr) & (lambda_r <= gamma_ra) & ...
         ~(test_el | test_s);
Psi(test_r) = 0.5 * lame(test_r) .* trace_E(test_r).^2 + ...
              shear(test_r) .* (eig_1(test_r).^2 + 0.5 * (eig_2(test_r) + eig_3(test_r)).^2) - ...
              0.5 * denom_r(test_r) .* (lambda_r(test_r).^2);  
              
% Return mapping: apex of the yield surface.
test_a = ~(test_el | test_s | test_l | test_r);
Psi(test_a) = 0.5 * bulk(test_a) .* trace_E(test_a).^2 - ...
              0.5 * denom_a(test_a) .* (lambda_a(test_a).^2);
          
end
