function [S,DS]=constitutive_problem_3D(E,c_bar,sin_phi,shear,bulk,lame)
%--------------------------------------------------------------------------
% constitutive_problem_3D computes the stress tensor and consistent tangent
% operator for a 3D elastic-perfectly plastic model based on the Mohr-Coulomb
% yield criterion with an associated flow rule.
%
% This function constructs both the constitutive operator (stress tensor S)
% and its consistent tangent derivative DS for each integration point, using a
% simplified static version of the elastic-plastic model. Different return
% mappings (elastic, smooth yield surface, left/right edge, apex) are determined
% by evaluating invariants and ordered eigenvalues of the trial strain tensor.
%
% INPUT ARGUMENTS:
%   E       - Current strain tensor at integration points, size (6, n_int),
%             where the strain components are ordered as [E11; E22; E33; E12; E13; E23].
%   c_bar   - Reduced cohesion parameter at integration points, size (1, n_int)
%   sin_phi - Sine of the reduced friction angle at integration points, size (1, n_int)
%   shear   - Shear modulus at integration points, size (1, n_int)
%   bulk    - Bulk modulus at integration points, size (1, n_int)
%   lame    - Lame's coefficient at integration points, size (1, n_int)
%
% OUTPUT:
%   S       - Stress tensor at integration points, size (6, n_int)
%   DS      - Consistent tangent operator at integration points, size (36, n_int)
%             (Only computed if nargout > 1)
%
%--------------------------------------------------------------------------

sin_psi=sin_phi;

%
% Linear constitutive operators
%
IDENT = diag([1, 1, 1, 1/2, 1/2, 1/2]) ; % identity operator
iota=[1;1;1;0;0;0];
VOL=iota*iota';
DEV=diag([1,1,1,1/2,1/2,1/2])-VOL/3;
ELAST=2*DEV(:)*shear+VOL(:)*bulk;   % size(ELAST)=(36,n_int)

%
% Trial strain and other auxiliary arrays:
%   E_tr   - trial strain tensors, size(E_tr)=(6,n_int)
%
n_int=size(E,2); % number of integration points
E_trial = E   ;  % strain representation
E_tr = IDENT*E_trial ;     % stress representation
E_square = ...             % square of the trial strain in stress repres.
    [ E_tr(1,:).^2         + E_tr(4,:).^2         + E_tr(6,:).^2
    E_tr(2,:).^2         + E_tr(4,:).^2         + E_tr(5,:).^2
    E_tr(3,:).^2         + E_tr(5,:).^2         + E_tr(6,:).^2
    E_tr(1,:).*E_tr(4,:) + E_tr(2,:).*E_tr(4,:) + E_tr(5,:).*E_tr(6,:)
    E_tr(4,:).*E_tr(6,:) + E_tr(2,:).*E_tr(5,:) + E_tr(3,:).*E_tr(5,:)
    E_tr(1,:).*E_tr(6,:) + E_tr(4,:).*E_tr(5,:) + E_tr(3,:).*E_tr(6,:) ];

%
% Invariants of the trial strain tensors (at integration points)
%
I1 = E_tr(1,:)+E_tr(2,:)+E_tr(3,:); % trace of E_trial
I2 = E_tr(1,:).*E_tr(2,:)+E_tr(1,:).*E_tr(3,:)+E_tr(2,:).*E_tr(3,:)- ...
    E_tr(4,:).^2 - E_tr(5,:).^2 - E_tr(6,:).^2;
I3 = E_tr(1,:).*E_tr(2,:).*E_tr(3,:) - E_tr(3,:).*E_tr(4,:).^2 - ...
    E_tr(2,:).*E_tr(6,:).^2 - E_tr(1,:).*E_tr(5,:).^2 + ...
    2*E_tr(4,:).*E_tr(5,:).*E_tr(6,:) ;

Q = max(0,(1/9)*((I1.^2) - 3*I2));
R = (1/54)*(-2*(I1.^3) + 9*(I1.*I2) - 27*I3);
test1 = (Q==0);
theta0 = zeros(1,n_int);
theta0(~test1) = R(~test1)./sqrt(Q(~test1).^3);
theta = acos(min(max(theta0,-1),1))/3; % Lode's angle

%
% Ordered eigenvalues of the trial strain tensors
%
eig_1 = -2*sqrt(Q).*cos(theta+2*pi/3) + I1/3 ;
eig_2 = -2*sqrt(Q).*cos(theta-2*pi/3) + I1/3 ;
eig_3 = -2*sqrt(Q).*cos(theta) + I1/3;

%
% Critical values defining decision criteria
%
f_tr=2*shear.*((1+sin_phi).*eig_1-(1-sin_phi).*eig_3)+ ...
    2*(lame.*sin_phi).*I1-c_bar;             % elast-plast interface
gamma_sl=(eig_1-eig_2)./(1+sin_psi);          % smooth-left interface
gamma_sr=(eig_2-eig_3)./(1-sin_psi);          % smooth-right interface
gamma_la=(eig_1+eig_2-2*eig_3)./(3-sin_psi);  % left-apex interface
gamma_ra=(2*eig_1-eig_2-eig_3)./(3+sin_psi);  % right-apex interface

%
% Candidates on plastic multipliers
%
denom_s=4*lame.*sin_phi.*sin_psi+4*shear.*(1+sin_phi.*sin_psi);
denom_l=4*lame.*sin_phi.*sin_psi+shear.*(1+sin_phi).*(1+sin_psi)+ ...
    2*shear.*(1-sin_phi).*(1-sin_psi);
denom_r=4*lame.*sin_phi.*sin_psi+2*shear.*(1+sin_phi).*(1+sin_psi)+ ...
    shear.*(1-sin_phi).*(1-sin_psi);
denom_a=4*bulk.*sin_phi.*sin_psi;  % denominators for each type of return

lambda_s=f_tr./denom_s ;
lambda_l=(shear.*((1+sin_phi).*(eig_1+eig_2)-2*(1-sin_phi).*eig_3)+ ...
    2*lame.*sin_phi.*I1-c_bar)./denom_l ;
lambda_r=(shear.*(2*(1+sin_phi).*eig_1-(1-sin_phi).*(eig_2+eig_3))+ ...
    2*lame.*sin_phi.*I1-c_bar)./denom_r ;
lambda_a=(2*bulk.*sin_phi.*I1-c_bar)./denom_a ;

%
% Determination of the stress array S
%

% initialization
S = zeros(6,n_int);

% elastic response
test_el = (f_tr<=0);
lame_el=lame(test_el);
shear_el=shear(test_el);
S(:,test_el)=  repmat(lame_el,6,1).*(VOL*E_trial(:,test_el))+...
    2*repmat(shear_el,6,1).*(IDENT*E_trial(:,test_el));

% return to the smooth portion of the yield surface
test_s = (lambda_s<=min(gamma_sl,gamma_sr))&(~test_el);
lame_s=lame(test_s);
shear_s=shear(test_s);
sin_phi_s=sin_phi(test_s);
sin_psi_s=sin_psi(test_s);
eig_1_s=eig_1(test_s);
eig_2_s=eig_2(test_s);
eig_3_s=eig_3(test_s);
I1_s=I1(test_s);
E_square_s=E_square(:,test_s);
E_tr_s=E_tr(:,test_s);
lambda_s=lambda_s(test_s);
% eigenprojections
denom_s1 = (eig_1_s-eig_2_s).*(eig_1_s-eig_3_s);
denom_s2 = (eig_2_s-eig_1_s).*(eig_2_s-eig_3_s);
denom_s3 = (eig_3_s-eig_1_s).*(eig_3_s-eig_2_s);
Eig_1_s = (ones(6,1)*(1./denom_s1)).* ( E_square_s -...
    (ones(6,1)*(eig_2_s+eig_3_s)).*E_tr_s +...
    iota*(eig_2_s.*eig_3_s) );
Eig_2_s = (ones(6,1)*(1./denom_s2)).* ( E_square_s -...
    (ones(6,1)*(eig_1_s+eig_3_s)).*E_tr_s +...
    iota*(eig_1_s.*eig_3_s) );
Eig_3_s = (ones(6,1)*(1./denom_s3)).* ( E_square_s -...
    (ones(6,1)*(eig_1_s+eig_2_s)).*E_tr_s +...
    iota*(eig_1_s.*eig_2_s) );
% principal stresses
sigma_1_s = lame_s.*I1_s+2*shear_s.*eig_1_s-...
    lambda_s.*(2*lame_s.*sin_psi_s+2*shear_s.*(1+sin_psi_s));
sigma_2_s = lame_s.*I1_s+2*shear_s.*eig_2_s-...
    lambda_s.*(2*lame_s.*sin_psi_s);
sigma_3_s = lame_s.*I1_s+2*shear_s.*eig_3_s-...
    lambda_s.*(2*lame_s.*sin_psi_s-2*shear_s.*(1-sin_psi_s));
% unknown stress tensors
S(:,test_s) = (ones(6,1)*sigma_1_s).*Eig_1_s+...
    (ones(6,1)*sigma_2_s).*Eig_2_s+...
    (ones(6,1)*sigma_3_s).*Eig_3_s;

% return to the left edge of the yield surface
test_l = (gamma_sl<gamma_sr)&(lambda_l>=gamma_sl)&...
    (lambda_l<=gamma_la)&(~(test_el|test_s));
lame_l=lame(test_l);
shear_l=shear(test_l);
sin_phi_l=sin_phi(test_l);
sin_psi_l=sin_psi(test_l);
nt_l=length(lame_l);
eig_1_l=eig_1(test_l);
eig_2_l=eig_2(test_l);
eig_3_l=eig_3(test_l);
I1_l=I1(test_l);
lambda_l=lambda_l(test_l);
E_square_l=E_square(:,test_l);
E_tr_l=E_tr(:,test_l);
% eigenprojections
denom_l3 = (eig_3_l-eig_1_l).*(eig_3_l-eig_2_l);
Eig_3_l  = (ones(6,1)*(1./denom_l3)).* ( E_square_l -...
    (ones(6,1)*(eig_1_l+eig_2_l)).*E_tr_l +...
    iota*(eig_1_l.*eig_2_l) );
Eig_12_l = [ones(3,nt_l); zeros(3,nt_l)] - Eig_3_l;
% principal stresses
sigma_1_l = lame_l.*I1_l+shear_l.*(eig_1_l+eig_2_l)-...
    lambda_l.*(2*lame_l.*sin_psi_l+shear_l.*(1+sin_psi_l));
sigma_3_l = lame_l.*I1_l+2*shear_l.*eig_3_l-...
    lambda_l.*(2*lame_l.*sin_psi_l-2*shear_l.*(1-sin_psi_l));
% unknown stress tensors
S(:,test_l) = (ones(6,1)*sigma_1_l).*Eig_12_l+...
    (ones(6,1)*sigma_3_l).*Eig_3_l;

% return to the right edge of the yield surface
test_r = (gamma_sl>gamma_sr)&(lambda_r>=gamma_sr)&...
    (lambda_r<=gamma_ra)&(~(test_el|test_s));
lame_r=lame(test_r);
shear_r=shear(test_r);
sin_phi_r=sin_phi(test_r);
sin_psi_r=sin_psi(test_r);
nt_r=length(lame_r);
eig_1_r=eig_1(test_r);
eig_2_r=eig_2(test_r);
eig_3_r=eig_3(test_r);
I1_r=I1(test_r);
lambda_r=lambda_r(test_r);
% eigenprojections
denom_r1 = (eig_1_r-eig_2_r).*(eig_1_r-eig_3_r);
Eig_1_r  = (ones(6,1)*(1./denom_r1)).* ( E_square(:,test_r) -...
    (ones(6,1)*(eig_2_r+eig_3_r)).*E_tr(:,test_r) +...
    iota*(eig_2_r.*eig_3_r) );
Eig_23_r = [ones(3,nt_r); zeros(3,nt_r)] - Eig_1_r ;
% principal stresses
sigma_1_r = lame_r.*I1_r+2*shear_r.*eig_1_r-...
    lambda_r.*(2*lame_r.*sin_psi_r+2*shear_r.*(1+sin_psi_r));
sigma_3_r = lame_r.*I1_r+shear_r.*(eig_2_r+eig_3_r)-...
    lambda_r.*(2*lame_r.*sin_psi_r-shear_r.*(1-sin_psi_r));
% unknown stress tensors
S(:,test_r) = (ones(6,1)*sigma_1_r).*Eig_1_r+...
    (ones(6,1)*sigma_3_r).*Eig_23_r;

% return to the apex of the yield surface
test_a=~(test_el|test_s|test_l|test_r);
lambda_a=lambda_a(test_a);
nt_a=length(lambda_a);
sigma_1_a = c_bar(test_a)./(2*sin_phi(test_a));
S(:,test_a) = iota*sigma_1_a;

%
% Determination of the tangential stiffness array DS
%

if nargout>1

    % initialization of the unknown array
    DS = zeros(36,n_int);

    % derivative of the square of the trial strain tensor(stress notation)
    DER_E_square = ... % 36*n_int array (fourth order tensors at int. points)
        [ 2*E_tr(1,:)   ; zeros(1,n_int); zeros(1,n_int); E_tr(4,:)                ; zeros(1,n_int)           ; E_tr(6,:)
        zeros(1,n_int); 2*E_tr(2,:)   ; zeros(1,n_int); E_tr(4,:)                ; E_tr(5,:)                ; zeros(1,n_int)
        zeros(1,n_int); zeros(1,n_int); 2*E_tr(3,:)   ; zeros(1,n_int)           ; E_tr(5,:)                ; E_tr(6,:)
        E_tr(4,:)     ; E_tr(4,:)     ; zeros(1,n_int); 0.5*(E_tr(1,:)+E_tr(2,:)); 0.5*E_tr(6,:)            ; 0.5*E_tr(5,:)
        zeros(1,n_int); E_tr(5,:)     ; E_tr(5,:)     ; 0.5*E_tr(6,:)            ; 0.5*(E_tr(2,:)+E_tr(3,:)); 0.5*E_tr(4,:)
        E_tr(6,:)     ; zeros(1,n_int); E_tr(6,:)     ; 0.5*E_tr(5,:)            ; 0.5*E_tr(4,:)            ; 0.5*(E_tr(1,:)+E_tr(3,:)) ];

    % elastic response
    DS(:,test_el)=ELAST(:,test_el);

    % return to the smooth portion of the yield surface
    % auxilliary arrays (36*n_int)
    E1_x_E1 = ...
        [ Eig_1_s(1,:).*Eig_1_s(1,:); Eig_1_s(2,:).*Eig_1_s(1,:); Eig_1_s(3,:).*Eig_1_s(1,:); Eig_1_s(4,:).*Eig_1_s(1,:); Eig_1_s(5,:).*Eig_1_s(1,:); Eig_1_s(6,:).*Eig_1_s(1,:)
        Eig_1_s(1,:).*Eig_1_s(2,:); Eig_1_s(2,:).*Eig_1_s(2,:); Eig_1_s(3,:).*Eig_1_s(2,:); Eig_1_s(4,:).*Eig_1_s(2,:); Eig_1_s(5,:).*Eig_1_s(2,:); Eig_1_s(6,:).*Eig_1_s(2,:)
        Eig_1_s(1,:).*Eig_1_s(3,:); Eig_1_s(2,:).*Eig_1_s(3,:); Eig_1_s(3,:).*Eig_1_s(3,:); Eig_1_s(4,:).*Eig_1_s(3,:); Eig_1_s(5,:).*Eig_1_s(3,:); Eig_1_s(6,:).*Eig_1_s(3,:)
        Eig_1_s(1,:).*Eig_1_s(4,:); Eig_1_s(2,:).*Eig_1_s(4,:); Eig_1_s(3,:).*Eig_1_s(4,:); Eig_1_s(4,:).*Eig_1_s(4,:); Eig_1_s(5,:).*Eig_1_s(4,:); Eig_1_s(6,:).*Eig_1_s(4,:)
        Eig_1_s(1,:).*Eig_1_s(5,:); Eig_1_s(2,:).*Eig_1_s(5,:); Eig_1_s(3,:).*Eig_1_s(5,:); Eig_1_s(4,:).*Eig_1_s(5,:); Eig_1_s(5,:).*Eig_1_s(5,:); Eig_1_s(6,:).*Eig_1_s(5,:)
        Eig_1_s(1,:).*Eig_1_s(6,:); Eig_1_s(2,:).*Eig_1_s(6,:); Eig_1_s(3,:).*Eig_1_s(6,:); Eig_1_s(4,:).*Eig_1_s(6,:); Eig_1_s(5,:).*Eig_1_s(6,:); Eig_1_s(6,:).*Eig_1_s(6,:) ];
    E2_x_E2 = ...
        [ Eig_2_s(1,:).*Eig_2_s(1,:); Eig_2_s(2,:).*Eig_2_s(1,:); Eig_2_s(3,:).*Eig_2_s(1,:); Eig_2_s(4,:).*Eig_2_s(1,:); Eig_2_s(5,:).*Eig_2_s(1,:); Eig_2_s(6,:).*Eig_2_s(1,:)
        Eig_2_s(1,:).*Eig_2_s(2,:); Eig_2_s(2,:).*Eig_2_s(2,:); Eig_2_s(3,:).*Eig_2_s(2,:); Eig_2_s(4,:).*Eig_2_s(2,:); Eig_2_s(5,:).*Eig_2_s(2,:); Eig_2_s(6,:).*Eig_2_s(2,:)
        Eig_2_s(1,:).*Eig_2_s(3,:); Eig_2_s(2,:).*Eig_2_s(3,:); Eig_2_s(3,:).*Eig_2_s(3,:); Eig_2_s(4,:).*Eig_2_s(3,:); Eig_2_s(5,:).*Eig_2_s(3,:); Eig_2_s(6,:).*Eig_2_s(3,:)
        Eig_2_s(1,:).*Eig_2_s(4,:); Eig_2_s(2,:).*Eig_2_s(4,:); Eig_2_s(3,:).*Eig_2_s(4,:); Eig_2_s(4,:).*Eig_2_s(4,:); Eig_2_s(5,:).*Eig_2_s(4,:); Eig_2_s(6,:).*Eig_2_s(4,:)
        Eig_2_s(1,:).*Eig_2_s(5,:); Eig_2_s(2,:).*Eig_2_s(5,:); Eig_2_s(3,:).*Eig_2_s(5,:); Eig_2_s(4,:).*Eig_2_s(5,:); Eig_2_s(5,:).*Eig_2_s(5,:); Eig_2_s(6,:).*Eig_2_s(5,:)
        Eig_2_s(1,:).*Eig_2_s(6,:); Eig_2_s(2,:).*Eig_2_s(6,:); Eig_2_s(3,:).*Eig_2_s(6,:); Eig_2_s(4,:).*Eig_2_s(6,:); Eig_2_s(5,:).*Eig_2_s(6,:); Eig_2_s(6,:).*Eig_2_s(6,:) ];
    E3_x_E3 = ...
        [ Eig_3_s(1,:).*Eig_3_s(1,:); Eig_3_s(2,:).*Eig_3_s(1,:); Eig_3_s(3,:).*Eig_3_s(1,:); Eig_3_s(4,:).*Eig_3_s(1,:); Eig_3_s(5,:).*Eig_3_s(1,:); Eig_3_s(6,:).*Eig_3_s(1,:)
        Eig_3_s(1,:).*Eig_3_s(2,:); Eig_3_s(2,:).*Eig_3_s(2,:); Eig_3_s(3,:).*Eig_3_s(2,:); Eig_3_s(4,:).*Eig_3_s(2,:); Eig_3_s(5,:).*Eig_3_s(2,:); Eig_3_s(6,:).*Eig_3_s(2,:)
        Eig_3_s(1,:).*Eig_3_s(3,:); Eig_3_s(2,:).*Eig_3_s(3,:); Eig_3_s(3,:).*Eig_3_s(3,:); Eig_3_s(4,:).*Eig_3_s(3,:); Eig_3_s(5,:).*Eig_3_s(3,:); Eig_3_s(6,:).*Eig_3_s(3,:)
        Eig_3_s(1,:).*Eig_3_s(4,:); Eig_3_s(2,:).*Eig_3_s(4,:); Eig_3_s(3,:).*Eig_3_s(4,:); Eig_3_s(4,:).*Eig_3_s(4,:); Eig_3_s(5,:).*Eig_3_s(4,:); Eig_3_s(6,:).*Eig_3_s(4,:)
        Eig_3_s(1,:).*Eig_3_s(5,:); Eig_3_s(2,:).*Eig_3_s(5,:); Eig_3_s(3,:).*Eig_3_s(5,:); Eig_3_s(4,:).*Eig_3_s(5,:); Eig_3_s(5,:).*Eig_3_s(5,:); Eig_3_s(6,:).*Eig_3_s(5,:)
        Eig_3_s(1,:).*Eig_3_s(6,:); Eig_3_s(2,:).*Eig_3_s(6,:); Eig_3_s(3,:).*Eig_3_s(6,:); Eig_3_s(4,:).*Eig_3_s(6,:); Eig_3_s(5,:).*Eig_3_s(6,:); Eig_3_s(6,:).*Eig_3_s(6,:) ];
    % derivatives of eigenprojections  (36*n_int)
    DER_E_square_s=DER_E_square(:,test_s);
    EIG_1_s = (ones(36,1)*(1./denom_s1)).*( ...
        DER_E_square_s - ...
        IDENT(:)*(eig_2_s+eig_3_s) - ...
        (ones(36,1)*(2*eig_1_s-eig_2_s-eig_3_s)).*E1_x_E1-...
        (ones(36,1)*(eig_2_s-eig_3_s)).*(E2_x_E2 - E3_x_E3) ...
        );
    EIG_2_s = (ones(36,1)*(1./denom_s2)).*( ...
        DER_E_square_s - ...
        IDENT(:)*(eig_1_s+eig_3_s) - ...
        (ones(36,1)*(2*eig_2_s-eig_1_s-eig_3_s)).*E2_x_E2-...
        (ones(36,1)*(eig_1_s-eig_3_s)).*(E1_x_E1 - E3_x_E3) ...
        );
    EIG_3_s = (ones(36,1)*(1./denom_s3)).*( ...
        DER_E_square_s - ...
        IDENT(:)*(eig_1_s+eig_2_s) - ...
        (ones(36,1)*(2*eig_3_s-eig_1_s-eig_2_s)).*E3_x_E3-...
        (ones(36,1)*(eig_1_s-eig_2_s)).*(E1_x_E1 - E2_x_E2) ...
        );
    % computation of the consistent tangent operators (36*n_int)
    Sder1_s=(ones(36,1)*sigma_1_s).*EIG_1_s+...
        (ones(36,1)*sigma_2_s).*EIG_2_s+...
        (ones(36,1)*sigma_3_s).*EIG_3_s;
    Sder2_s=VOL(:)*lame_s;
    Sder3_s=2*repmat(shear_s,36,1).*( E1_x_E1 + E2_x_E2 + E3_x_E3 ) ;
    D_phi_s = 2*repmat(shear_s,6,1).*(repmat(1+sin_phi_s,6,1).*Eig_1_s-repmat(1-sin_phi_s,6,1).*Eig_3_s)...
        +2*iota*(lame_s.*sin_phi_s);
    D_psi_s = 2*repmat(shear_s,6,1).*(repmat(1+sin_psi_s,6,1).*Eig_1_s-repmat(1-sin_psi_s,6,1).*Eig_3_s)...
        +2*iota*(lame_s.*sin_psi_s);
    Sder4_s = ...
        [ D_psi_s(1,:).*D_phi_s(1,:); D_psi_s(2,:).*D_phi_s(1,:); D_psi_s(3,:).*D_phi_s(1,:); D_psi_s(4,:).*D_phi_s(1,:); D_psi_s(5,:).*D_phi_s(1,:); D_psi_s(6,:).*D_phi_s(1,:)
        D_psi_s(1,:).*D_phi_s(2,:); D_psi_s(2,:).*D_phi_s(2,:); D_psi_s(3,:).*D_phi_s(2,:); D_psi_s(4,:).*D_phi_s(2,:); D_psi_s(5,:).*D_phi_s(2,:); D_psi_s(6,:).*D_phi_s(2,:)
        D_psi_s(1,:).*D_phi_s(3,:); D_psi_s(2,:).*D_phi_s(3,:); D_psi_s(3,:).*D_phi_s(3,:); D_psi_s(4,:).*D_phi_s(3,:); D_psi_s(5,:).*D_phi_s(3,:); D_psi_s(6,:).*D_phi_s(3,:)
        D_psi_s(1,:).*D_phi_s(4,:); D_psi_s(2,:).*D_phi_s(4,:); D_psi_s(3,:).*D_phi_s(4,:); D_psi_s(4,:).*D_phi_s(4,:); D_psi_s(5,:).*D_phi_s(4,:); D_psi_s(6,:).*D_phi_s(4,:)
        D_psi_s(1,:).*D_phi_s(5,:); D_psi_s(2,:).*D_phi_s(5,:); D_psi_s(3,:).*D_phi_s(5,:); D_psi_s(4,:).*D_phi_s(5,:); D_psi_s(5,:).*D_phi_s(5,:); D_psi_s(6,:).*D_phi_s(5,:)
        D_psi_s(1,:).*D_phi_s(6,:); D_psi_s(2,:).*D_phi_s(6,:); D_psi_s(3,:).*D_phi_s(6,:); D_psi_s(4,:).*D_phi_s(6,:); D_psi_s(5,:).*D_phi_s(6,:); D_psi_s(6,:).*D_phi_s(6,:) ]./repmat(denom_s(test_s),36,1);
    DS(:,test_s)=Sder1_s+Sder2_s+Sder3_s-Sder4_s;

    % return to the left edge of the yield surface
    % auxilliary arrays (36*n_int)
    E3_x_E3 = ...
        [ Eig_3_l(1,:).*Eig_3_l(1,:); Eig_3_l(2,:).*Eig_3_l(1,:); Eig_3_l(3,:).*Eig_3_l(1,:); Eig_3_l(4,:).*Eig_3_l(1,:); Eig_3_l(5,:).*Eig_3_l(1,:); Eig_3_l(6,:).*Eig_3_l(1,:)
        Eig_3_l(1,:).*Eig_3_l(2,:); Eig_3_l(2,:).*Eig_3_l(2,:); Eig_3_l(3,:).*Eig_3_l(2,:); Eig_3_l(4,:).*Eig_3_l(2,:); Eig_3_l(5,:).*Eig_3_l(2,:); Eig_3_l(6,:).*Eig_3_l(2,:)
        Eig_3_l(1,:).*Eig_3_l(3,:); Eig_3_l(2,:).*Eig_3_l(3,:); Eig_3_l(3,:).*Eig_3_l(3,:); Eig_3_l(4,:).*Eig_3_l(3,:); Eig_3_l(5,:).*Eig_3_l(3,:); Eig_3_l(6,:).*Eig_3_l(3,:)
        Eig_3_l(1,:).*Eig_3_l(4,:); Eig_3_l(2,:).*Eig_3_l(4,:); Eig_3_l(3,:).*Eig_3_l(4,:); Eig_3_l(4,:).*Eig_3_l(4,:); Eig_3_l(5,:).*Eig_3_l(4,:); Eig_3_l(6,:).*Eig_3_l(4,:)
        Eig_3_l(1,:).*Eig_3_l(5,:); Eig_3_l(2,:).*Eig_3_l(5,:); Eig_3_l(3,:).*Eig_3_l(5,:); Eig_3_l(4,:).*Eig_3_l(5,:); Eig_3_l(5,:).*Eig_3_l(5,:); Eig_3_l(6,:).*Eig_3_l(5,:)
        Eig_3_l(1,:).*Eig_3_l(6,:); Eig_3_l(2,:).*Eig_3_l(6,:); Eig_3_l(3,:).*Eig_3_l(6,:); Eig_3_l(4,:).*Eig_3_l(6,:); Eig_3_l(5,:).*Eig_3_l(6,:); Eig_3_l(6,:).*Eig_3_l(6,:) ];
    E12_x_E12 = ...
        [ Eig_12_l(1,:).*Eig_12_l(1,:); Eig_12_l(2,:).*Eig_12_l(1,:); Eig_12_l(3,:).*Eig_12_l(1,:); Eig_12_l(4,:).*Eig_12_l(1,:); Eig_12_l(5,:).*Eig_12_l(1,:); Eig_12_l(6,:).*Eig_12_l(1,:)
        Eig_12_l(1,:).*Eig_12_l(2,:); Eig_12_l(2,:).*Eig_12_l(2,:); Eig_12_l(3,:).*Eig_12_l(2,:); Eig_12_l(4,:).*Eig_12_l(2,:); Eig_12_l(5,:).*Eig_12_l(2,:); Eig_12_l(6,:).*Eig_12_l(2,:)
        Eig_12_l(1,:).*Eig_12_l(3,:); Eig_12_l(2,:).*Eig_12_l(3,:); Eig_12_l(3,:).*Eig_12_l(3,:); Eig_12_l(4,:).*Eig_12_l(3,:); Eig_12_l(5,:).*Eig_12_l(3,:); Eig_12_l(6,:).*Eig_12_l(3,:)
        Eig_12_l(1,:).*Eig_12_l(4,:); Eig_12_l(2,:).*Eig_12_l(4,:); Eig_12_l(3,:).*Eig_12_l(4,:); Eig_12_l(4,:).*Eig_12_l(4,:); Eig_12_l(5,:).*Eig_12_l(4,:); Eig_12_l(6,:).*Eig_12_l(4,:)
        Eig_12_l(1,:).*Eig_12_l(5,:); Eig_12_l(2,:).*Eig_12_l(5,:); Eig_12_l(3,:).*Eig_12_l(5,:); Eig_12_l(4,:).*Eig_12_l(5,:); Eig_12_l(5,:).*Eig_12_l(5,:); Eig_12_l(6,:).*Eig_12_l(5,:)
        Eig_12_l(1,:).*Eig_12_l(6,:); Eig_12_l(2,:).*Eig_12_l(6,:); Eig_12_l(3,:).*Eig_12_l(6,:); Eig_12_l(4,:).*Eig_12_l(6,:); Eig_12_l(5,:).*Eig_12_l(6,:); Eig_12_l(6,:).*Eig_12_l(6,:) ];
    E_tr_l = E_tr(:,test_l) ;
    E12_x_Etr = ...
        [ Eig_12_l(1,:).*E_tr_l(1,:); Eig_12_l(2,:).*E_tr_l(1,:); Eig_12_l(3,:).*E_tr_l(1,:); Eig_12_l(4,:).*E_tr_l(1,:); Eig_12_l(5,:).*E_tr_l(1,:); Eig_12_l(6,:).*E_tr_l(1,:)
        Eig_12_l(1,:).*E_tr_l(2,:); Eig_12_l(2,:).*E_tr_l(2,:); Eig_12_l(3,:).*E_tr_l(2,:); Eig_12_l(4,:).*E_tr_l(2,:); Eig_12_l(5,:).*E_tr_l(2,:); Eig_12_l(6,:).*E_tr_l(2,:)
        Eig_12_l(1,:).*E_tr_l(3,:); Eig_12_l(2,:).*E_tr_l(3,:); Eig_12_l(3,:).*E_tr_l(3,:); Eig_12_l(4,:).*E_tr_l(3,:); Eig_12_l(5,:).*E_tr_l(3,:); Eig_12_l(6,:).*E_tr_l(3,:)
        Eig_12_l(1,:).*E_tr_l(4,:); Eig_12_l(2,:).*E_tr_l(4,:); Eig_12_l(3,:).*E_tr_l(4,:); Eig_12_l(4,:).*E_tr_l(4,:); Eig_12_l(5,:).*E_tr_l(4,:); Eig_12_l(6,:).*E_tr_l(4,:)
        Eig_12_l(1,:).*E_tr_l(5,:); Eig_12_l(2,:).*E_tr_l(5,:); Eig_12_l(3,:).*E_tr_l(5,:); Eig_12_l(4,:).*E_tr_l(5,:); Eig_12_l(5,:).*E_tr_l(5,:); Eig_12_l(6,:).*E_tr_l(5,:)
        Eig_12_l(1,:).*E_tr_l(6,:); Eig_12_l(2,:).*E_tr_l(6,:); Eig_12_l(3,:).*E_tr_l(6,:); Eig_12_l(4,:).*E_tr_l(6,:); Eig_12_l(5,:).*E_tr_l(6,:); Eig_12_l(6,:).*E_tr_l(6,:) ];
    Etr_x_E12 = ...
        [ E_tr_l(1,:).*Eig_12_l(1,:); E_tr_l(2,:).*Eig_12_l(1,:); E_tr_l(3,:).*Eig_12_l(1,:); E_tr_l(4,:).*Eig_12_l(1,:); E_tr_l(5,:).*Eig_12_l(1,:); E_tr_l(6,:).*Eig_12_l(1,:)
        E_tr_l(1,:).*Eig_12_l(2,:); E_tr_l(2,:).*Eig_12_l(2,:); E_tr_l(3,:).*Eig_12_l(2,:); E_tr_l(4,:).*Eig_12_l(2,:); E_tr_l(5,:).*Eig_12_l(2,:); E_tr_l(6,:).*Eig_12_l(2,:)
        E_tr_l(1,:).*Eig_12_l(3,:); E_tr_l(2,:).*Eig_12_l(3,:); E_tr_l(3,:).*Eig_12_l(3,:); E_tr_l(4,:).*Eig_12_l(3,:); E_tr_l(5,:).*Eig_12_l(3,:); E_tr_l(6,:).*Eig_12_l(3,:)
        E_tr_l(1,:).*Eig_12_l(4,:); E_tr_l(2,:).*Eig_12_l(4,:); E_tr_l(3,:).*Eig_12_l(4,:); E_tr_l(4,:).*Eig_12_l(4,:); E_tr_l(5,:).*Eig_12_l(4,:); E_tr_l(6,:).*Eig_12_l(4,:)
        E_tr_l(1,:).*Eig_12_l(5,:); E_tr_l(2,:).*Eig_12_l(5,:); E_tr_l(3,:).*Eig_12_l(5,:); E_tr_l(4,:).*Eig_12_l(5,:); E_tr_l(5,:).*Eig_12_l(5,:); E_tr_l(6,:).*Eig_12_l(5,:)
        E_tr_l(1,:).*Eig_12_l(6,:); E_tr_l(2,:).*Eig_12_l(6,:); E_tr_l(3,:).*Eig_12_l(6,:); E_tr_l(4,:).*Eig_12_l(6,:); E_tr_l(5,:).*Eig_12_l(6,:); E_tr_l(6,:).*Eig_12_l(6,:) ];
    E12_x_E3 = ...
        [ Eig_12_l(1,:).*Eig_3_l(1,:); Eig_12_l(2,:).*Eig_3_l(1,:); Eig_12_l(3,:).*Eig_3_l(1,:); Eig_12_l(4,:).*Eig_3_l(1,:); Eig_12_l(5,:).*Eig_3_l(1,:); Eig_12_l(6,:).*Eig_3_l(1,:)
        Eig_12_l(1,:).*Eig_3_l(2,:); Eig_12_l(2,:).*Eig_3_l(2,:); Eig_12_l(3,:).*Eig_3_l(2,:); Eig_12_l(4,:).*Eig_3_l(2,:); Eig_12_l(5,:).*Eig_3_l(2,:); Eig_12_l(6,:).*Eig_3_l(2,:)
        Eig_12_l(1,:).*Eig_3_l(3,:); Eig_12_l(2,:).*Eig_3_l(3,:); Eig_12_l(3,:).*Eig_3_l(3,:); Eig_12_l(4,:).*Eig_3_l(3,:); Eig_12_l(5,:).*Eig_3_l(3,:); Eig_12_l(6,:).*Eig_3_l(3,:)
        Eig_12_l(1,:).*Eig_3_l(4,:); Eig_12_l(2,:).*Eig_3_l(4,:); Eig_12_l(3,:).*Eig_3_l(4,:); Eig_12_l(4,:).*Eig_3_l(4,:); Eig_12_l(5,:).*Eig_3_l(4,:); Eig_12_l(6,:).*Eig_3_l(4,:)
        Eig_12_l(1,:).*Eig_3_l(5,:); Eig_12_l(2,:).*Eig_3_l(5,:); Eig_12_l(3,:).*Eig_3_l(5,:); Eig_12_l(4,:).*Eig_3_l(5,:); Eig_12_l(5,:).*Eig_3_l(5,:); Eig_12_l(6,:).*Eig_3_l(5,:)
        Eig_12_l(1,:).*Eig_3_l(6,:); Eig_12_l(2,:).*Eig_3_l(6,:); Eig_12_l(3,:).*Eig_3_l(6,:); Eig_12_l(4,:).*Eig_3_l(6,:); Eig_12_l(5,:).*Eig_3_l(6,:); Eig_12_l(6,:).*Eig_3_l(6,:) ];
    E3_x_E12 = ...
        [ Eig_3_l(1,:).*Eig_12_l(1,:); Eig_3_l(2,:).*Eig_12_l(1,:); Eig_3_l(3,:).*Eig_12_l(1,:); Eig_3_l(4,:).*Eig_12_l(1,:); Eig_3_l(5,:).*Eig_12_l(1,:); Eig_3_l(6,:).*Eig_12_l(1,:)
        Eig_3_l(1,:).*Eig_12_l(2,:); Eig_3_l(2,:).*Eig_12_l(2,:); Eig_3_l(3,:).*Eig_12_l(2,:); Eig_3_l(4,:).*Eig_12_l(2,:); Eig_3_l(5,:).*Eig_12_l(2,:); Eig_3_l(6,:).*Eig_12_l(2,:)
        Eig_3_l(1,:).*Eig_12_l(3,:); Eig_3_l(2,:).*Eig_12_l(3,:); Eig_3_l(3,:).*Eig_12_l(3,:); Eig_3_l(4,:).*Eig_12_l(3,:); Eig_3_l(5,:).*Eig_12_l(3,:); Eig_3_l(6,:).*Eig_12_l(3,:)
        Eig_3_l(1,:).*Eig_12_l(4,:); Eig_3_l(2,:).*Eig_12_l(4,:); Eig_3_l(3,:).*Eig_12_l(4,:); Eig_3_l(4,:).*Eig_12_l(4,:); Eig_3_l(5,:).*Eig_12_l(4,:); Eig_3_l(6,:).*Eig_12_l(4,:)
        Eig_3_l(1,:).*Eig_12_l(5,:); Eig_3_l(2,:).*Eig_12_l(5,:); Eig_3_l(3,:).*Eig_12_l(5,:); Eig_3_l(4,:).*Eig_12_l(5,:); Eig_3_l(5,:).*Eig_12_l(5,:); Eig_3_l(6,:).*Eig_12_l(5,:)
        Eig_3_l(1,:).*Eig_12_l(6,:); Eig_3_l(2,:).*Eig_12_l(6,:); Eig_3_l(3,:).*Eig_12_l(6,:); Eig_3_l(4,:).*Eig_12_l(6,:); Eig_3_l(5,:).*Eig_12_l(6,:); Eig_3_l(6,:).*Eig_12_l(6,:) ];
    % derivative of the third eigenprojections (36*n_int)
    EIG_3_l = (ones(36,1)*(1./denom_l3)).*( ...
        DER_E_square(:,test_l)- ...
        IDENT(:)*(eig_1_l+eig_2_l) - ...
        (Etr_x_E12 + E12_x_Etr) + ...
        (ones(36,1)*(eig_1_l+eig_2_l)).*E12_x_E12 + ...
        (ones(36,1)*(eig_1_l+eig_2_l-2*eig_3_l)).*E3_x_E3 +...
        (ones(36,1)*(eig_3_l)).*(E12_x_E3 + E3_x_E12) ...
        );
    % computation of the consistent tangent operators (36*n_int)
    Sder1_l=(ones(36,1)*(sigma_3_l-sigma_1_l)).*EIG_3_l;
    Sder2_l=VOL(:)*lame_l;
    Sder3_l=repmat(shear_l,36,1).*( E12_x_E12 + 2*E3_x_E3 );
    D_phi_l=repmat(shear_l,6,1).*(repmat(1+sin_phi_l,6,1).*Eig_12_l-2*repmat(1-sin_phi_l,6,1).*Eig_3_l)...
        +2*iota*(lame_l.*sin_phi_l);
    D_psi_l=repmat(shear_l,6,1).*(repmat(1+sin_psi_l,6,1).*Eig_12_l-2*repmat(1-sin_psi_l,6,1).*Eig_3_l)...
        +2*iota*(lame_l.*sin_psi_l);
    Sder4_l = ...
        [ D_psi_l(1,:).*D_phi_l(1,:); D_psi_l(2,:).*D_phi_l(1,:); D_psi_l(3,:).*D_phi_l(1,:); D_psi_l(4,:).*D_phi_l(1,:); D_psi_l(5,:).*D_phi_l(1,:); D_psi_l(6,:).*D_phi_l(1,:)
        D_psi_l(1,:).*D_phi_l(2,:); D_psi_l(2,:).*D_phi_l(2,:); D_psi_l(3,:).*D_phi_l(2,:); D_psi_l(4,:).*D_phi_l(2,:); D_psi_l(5,:).*D_phi_l(2,:); D_psi_l(6,:).*D_phi_l(2,:)
        D_psi_l(1,:).*D_phi_l(3,:); D_psi_l(2,:).*D_phi_l(3,:); D_psi_l(3,:).*D_phi_l(3,:); D_psi_l(4,:).*D_phi_l(3,:); D_psi_l(5,:).*D_phi_l(3,:); D_psi_l(6,:).*D_phi_l(3,:)
        D_psi_l(1,:).*D_phi_l(4,:); D_psi_l(2,:).*D_phi_l(4,:); D_psi_l(3,:).*D_phi_l(4,:); D_psi_l(4,:).*D_phi_l(4,:); D_psi_l(5,:).*D_phi_l(4,:); D_psi_l(6,:).*D_phi_l(4,:)
        D_psi_l(1,:).*D_phi_l(5,:); D_psi_l(2,:).*D_phi_l(5,:); D_psi_l(3,:).*D_phi_l(5,:); D_psi_l(4,:).*D_phi_l(5,:); D_psi_l(5,:).*D_phi_l(5,:); D_psi_l(6,:).*D_phi_l(5,:)
        D_psi_l(1,:).*D_phi_l(6,:); D_psi_l(2,:).*D_phi_l(6,:); D_psi_l(3,:).*D_phi_l(6,:); D_psi_l(4,:).*D_phi_l(6,:); D_psi_l(5,:).*D_phi_l(6,:); D_psi_l(6,:).*D_phi_l(6,:) ]./repmat(denom_l(test_l),36,1);
    DS(:,test_l)=Sder1_l+Sder2_l+Sder3_l-Sder4_l;

    % return to the right edge of the yield surface
    % auxilliary arrays (36*n_int)
    E1_x_E1 = ...
        [ Eig_1_r(1,:).*Eig_1_r(1,:); Eig_1_r(2,:).*Eig_1_r(1,:); Eig_1_r(3,:).*Eig_1_r(1,:); Eig_1_r(4,:).*Eig_1_r(1,:); Eig_1_r(5,:).*Eig_1_r(1,:); Eig_1_r(6,:).*Eig_1_r(1,:)
        Eig_1_r(1,:).*Eig_1_r(2,:); Eig_1_r(2,:).*Eig_1_r(2,:); Eig_1_r(3,:).*Eig_1_r(2,:); Eig_1_r(4,:).*Eig_1_r(2,:); Eig_1_r(5,:).*Eig_1_r(2,:); Eig_1_r(6,:).*Eig_1_r(2,:)
        Eig_1_r(1,:).*Eig_1_r(3,:); Eig_1_r(2,:).*Eig_1_r(3,:); Eig_1_r(3,:).*Eig_1_r(3,:); Eig_1_r(4,:).*Eig_1_r(3,:); Eig_1_r(5,:).*Eig_1_r(3,:); Eig_1_r(6,:).*Eig_1_r(3,:)
        Eig_1_r(1,:).*Eig_1_r(4,:); Eig_1_r(2,:).*Eig_1_r(4,:); Eig_1_r(3,:).*Eig_1_r(4,:); Eig_1_r(4,:).*Eig_1_r(4,:); Eig_1_r(5,:).*Eig_1_r(4,:); Eig_1_r(6,:).*Eig_1_r(4,:)
        Eig_1_r(1,:).*Eig_1_r(5,:); Eig_1_r(2,:).*Eig_1_r(5,:); Eig_1_r(3,:).*Eig_1_r(5,:); Eig_1_r(4,:).*Eig_1_r(5,:); Eig_1_r(5,:).*Eig_1_r(5,:); Eig_1_r(6,:).*Eig_1_r(5,:)
        Eig_1_r(1,:).*Eig_1_r(6,:); Eig_1_r(2,:).*Eig_1_r(6,:); Eig_1_r(3,:).*Eig_1_r(6,:); Eig_1_r(4,:).*Eig_1_r(6,:); Eig_1_r(5,:).*Eig_1_r(6,:); Eig_1_r(6,:).*Eig_1_r(6,:) ];
    E23_x_E23 = ...
        [ Eig_23_r(1,:).*Eig_23_r(1,:); Eig_23_r(2,:).*Eig_23_r(1,:); Eig_23_r(3,:).*Eig_23_r(1,:); Eig_23_r(4,:).*Eig_23_r(1,:); Eig_23_r(5,:).*Eig_23_r(1,:); Eig_23_r(6,:).*Eig_23_r(1,:)
        Eig_23_r(1,:).*Eig_23_r(2,:); Eig_23_r(2,:).*Eig_23_r(2,:); Eig_23_r(3,:).*Eig_23_r(2,:); Eig_23_r(4,:).*Eig_23_r(2,:); Eig_23_r(5,:).*Eig_23_r(2,:); Eig_23_r(6,:).*Eig_23_r(2,:)
        Eig_23_r(1,:).*Eig_23_r(3,:); Eig_23_r(2,:).*Eig_23_r(3,:); Eig_23_r(3,:).*Eig_23_r(3,:); Eig_23_r(4,:).*Eig_23_r(3,:); Eig_23_r(5,:).*Eig_23_r(3,:); Eig_23_r(6,:).*Eig_23_r(3,:)
        Eig_23_r(1,:).*Eig_23_r(4,:); Eig_23_r(2,:).*Eig_23_r(4,:); Eig_23_r(3,:).*Eig_23_r(4,:); Eig_23_r(4,:).*Eig_23_r(4,:); Eig_23_r(5,:).*Eig_23_r(4,:); Eig_23_r(6,:).*Eig_23_r(4,:)
        Eig_23_r(1,:).*Eig_23_r(5,:); Eig_23_r(2,:).*Eig_23_r(5,:); Eig_23_r(3,:).*Eig_23_r(5,:); Eig_23_r(4,:).*Eig_23_r(5,:); Eig_23_r(5,:).*Eig_23_r(5,:); Eig_23_r(6,:).*Eig_23_r(5,:)
        Eig_23_r(1,:).*Eig_23_r(6,:); Eig_23_r(2,:).*Eig_23_r(6,:); Eig_23_r(3,:).*Eig_23_r(6,:); Eig_23_r(4,:).*Eig_23_r(6,:); Eig_23_r(5,:).*Eig_23_r(6,:); Eig_23_r(6,:).*Eig_23_r(6,:) ];
    E_tr_r = E_tr(:,test_r) ;
    E23_x_Etr = ...
        [ Eig_23_r(1,:).*E_tr_r(1,:); Eig_23_r(2,:).*E_tr_r(1,:); Eig_23_r(3,:).*E_tr_r(1,:); Eig_23_r(4,:).*E_tr_r(1,:); Eig_23_r(5,:).*E_tr_r(1,:); Eig_23_r(6,:).*E_tr_r(1,:)
        Eig_23_r(1,:).*E_tr_r(2,:); Eig_23_r(2,:).*E_tr_r(2,:); Eig_23_r(3,:).*E_tr_r(2,:); Eig_23_r(4,:).*E_tr_r(2,:); Eig_23_r(5,:).*E_tr_r(2,:); Eig_23_r(6,:).*E_tr_r(2,:)
        Eig_23_r(1,:).*E_tr_r(3,:); Eig_23_r(2,:).*E_tr_r(3,:); Eig_23_r(3,:).*E_tr_r(3,:); Eig_23_r(4,:).*E_tr_r(3,:); Eig_23_r(5,:).*E_tr_r(3,:); Eig_23_r(6,:).*E_tr_r(3,:)
        Eig_23_r(1,:).*E_tr_r(4,:); Eig_23_r(2,:).*E_tr_r(4,:); Eig_23_r(3,:).*E_tr_r(4,:); Eig_23_r(4,:).*E_tr_r(4,:); Eig_23_r(5,:).*E_tr_r(4,:); Eig_23_r(6,:).*E_tr_r(4,:)
        Eig_23_r(1,:).*E_tr_r(5,:); Eig_23_r(2,:).*E_tr_r(5,:); Eig_23_r(3,:).*E_tr_r(5,:); Eig_23_r(4,:).*E_tr_r(5,:); Eig_23_r(5,:).*E_tr_r(5,:); Eig_23_r(6,:).*E_tr_r(5,:)
        Eig_23_r(1,:).*E_tr_r(6,:); Eig_23_r(2,:).*E_tr_r(6,:); Eig_23_r(3,:).*E_tr_r(6,:); Eig_23_r(4,:).*E_tr_r(6,:); Eig_23_r(5,:).*E_tr_r(6,:); Eig_23_r(6,:).*E_tr_r(6,:) ];
    Etr_x_E23 = ...
        [ E_tr_r(1,:).*Eig_23_r(1,:); E_tr_r(2,:).*Eig_23_r(1,:); E_tr_r(3,:).*Eig_23_r(1,:); E_tr_r(4,:).*Eig_23_r(1,:); E_tr_r(5,:).*Eig_23_r(1,:); E_tr_r(6,:).*Eig_23_r(1,:)
        E_tr_r(1,:).*Eig_23_r(2,:); E_tr_r(2,:).*Eig_23_r(2,:); E_tr_r(3,:).*Eig_23_r(2,:); E_tr_r(4,:).*Eig_23_r(2,:); E_tr_r(5,:).*Eig_23_r(2,:); E_tr_r(6,:).*Eig_23_r(2,:)
        E_tr_r(1,:).*Eig_23_r(3,:); E_tr_r(2,:).*Eig_23_r(3,:); E_tr_r(3,:).*Eig_23_r(3,:); E_tr_r(4,:).*Eig_23_r(3,:); E_tr_r(5,:).*Eig_23_r(3,:); E_tr_r(6,:).*Eig_23_r(3,:)
        E_tr_r(1,:).*Eig_23_r(4,:); E_tr_r(2,:).*Eig_23_r(4,:); E_tr_r(3,:).*Eig_23_r(4,:); E_tr_r(4,:).*Eig_23_r(4,:); E_tr_r(5,:).*Eig_23_r(4,:); E_tr_r(6,:).*Eig_23_r(4,:)
        E_tr_r(1,:).*Eig_23_r(5,:); E_tr_r(2,:).*Eig_23_r(5,:); E_tr_r(3,:).*Eig_23_r(5,:); E_tr_r(4,:).*Eig_23_r(5,:); E_tr_r(5,:).*Eig_23_r(5,:); E_tr_r(6,:).*Eig_23_r(5,:)
        E_tr_r(1,:).*Eig_23_r(6,:); E_tr_r(2,:).*Eig_23_r(6,:); E_tr_r(3,:).*Eig_23_r(6,:); E_tr_r(4,:).*Eig_23_r(6,:); E_tr_r(5,:).*Eig_23_r(6,:); E_tr_r(6,:).*Eig_23_r(6,:) ];
    E23_x_E1 = ...
        [ Eig_23_r(1,:).*Eig_1_r(1,:); Eig_23_r(2,:).*Eig_1_r(1,:); Eig_23_r(3,:).*Eig_1_r(1,:); Eig_23_r(4,:).*Eig_1_r(1,:); Eig_23_r(5,:).*Eig_1_r(1,:); Eig_23_r(6,:).*Eig_1_r(1,:)
        Eig_23_r(1,:).*Eig_1_r(2,:); Eig_23_r(2,:).*Eig_1_r(2,:); Eig_23_r(3,:).*Eig_1_r(2,:); Eig_23_r(4,:).*Eig_1_r(2,:); Eig_23_r(5,:).*Eig_1_r(2,:); Eig_23_r(6,:).*Eig_1_r(2,:)
        Eig_23_r(1,:).*Eig_1_r(3,:); Eig_23_r(2,:).*Eig_1_r(3,:); Eig_23_r(3,:).*Eig_1_r(3,:); Eig_23_r(4,:).*Eig_1_r(3,:); Eig_23_r(5,:).*Eig_1_r(3,:); Eig_23_r(6,:).*Eig_1_r(3,:)
        Eig_23_r(1,:).*Eig_1_r(4,:); Eig_23_r(2,:).*Eig_1_r(4,:); Eig_23_r(3,:).*Eig_1_r(4,:); Eig_23_r(4,:).*Eig_1_r(4,:); Eig_23_r(5,:).*Eig_1_r(4,:); Eig_23_r(6,:).*Eig_1_r(4,:)
        Eig_23_r(1,:).*Eig_1_r(5,:); Eig_23_r(2,:).*Eig_1_r(5,:); Eig_23_r(3,:).*Eig_1_r(5,:); Eig_23_r(4,:).*Eig_1_r(5,:); Eig_23_r(5,:).*Eig_1_r(5,:); Eig_23_r(6,:).*Eig_1_r(5,:)
        Eig_23_r(1,:).*Eig_1_r(6,:); Eig_23_r(2,:).*Eig_1_r(6,:); Eig_23_r(3,:).*Eig_1_r(6,:); Eig_23_r(4,:).*Eig_1_r(6,:); Eig_23_r(5,:).*Eig_1_r(6,:); Eig_23_r(6,:).*Eig_1_r(6,:) ];
    E1_x_E23 = ...
        [ Eig_1_r(1,:).*Eig_23_r(1,:); Eig_1_r(2,:).*Eig_23_r(1,:); Eig_1_r(3,:).*Eig_23_r(1,:); Eig_1_r(4,:).*Eig_23_r(1,:); Eig_1_r(5,:).*Eig_23_r(1,:); Eig_1_r(6,:).*Eig_23_r(1,:)
        Eig_1_r(1,:).*Eig_23_r(2,:); Eig_1_r(2,:).*Eig_23_r(2,:); Eig_1_r(3,:).*Eig_23_r(2,:); Eig_1_r(4,:).*Eig_23_r(2,:); Eig_1_r(5,:).*Eig_23_r(2,:); Eig_1_r(6,:).*Eig_23_r(2,:)
        Eig_1_r(1,:).*Eig_23_r(3,:); Eig_1_r(2,:).*Eig_23_r(3,:); Eig_1_r(3,:).*Eig_23_r(3,:); Eig_1_r(4,:).*Eig_23_r(3,:); Eig_1_r(5,:).*Eig_23_r(3,:); Eig_1_r(6,:).*Eig_23_r(3,:)
        Eig_1_r(1,:).*Eig_23_r(4,:); Eig_1_r(2,:).*Eig_23_r(4,:); Eig_1_r(3,:).*Eig_23_r(4,:); Eig_1_r(4,:).*Eig_23_r(4,:); Eig_1_r(5,:).*Eig_23_r(4,:); Eig_1_r(6,:).*Eig_23_r(4,:)
        Eig_1_r(1,:).*Eig_23_r(5,:); Eig_1_r(2,:).*Eig_23_r(5,:); Eig_1_r(3,:).*Eig_23_r(5,:); Eig_1_r(4,:).*Eig_23_r(5,:); Eig_1_r(5,:).*Eig_23_r(5,:); Eig_1_r(6,:).*Eig_23_r(5,:)
        Eig_1_r(1,:).*Eig_23_r(6,:); Eig_1_r(2,:).*Eig_23_r(6,:); Eig_1_r(3,:).*Eig_23_r(6,:); Eig_1_r(4,:).*Eig_23_r(6,:); Eig_1_r(5,:).*Eig_23_r(6,:); Eig_1_r(6,:).*Eig_23_r(6,:) ];
    % derivative of the first eigenprojections (36*n_int)
    EIG_1_r = (ones(36,1)*(1./denom_r1)).*( DER_E_square(:,test_r) - ...
        IDENT(:)*(eig_2_r+eig_3_r) - ...
        Etr_x_E23 - E23_x_Etr + (ones(36,1)*(eig_2_r+eig_3_r)).*E23_x_E23 + ...
        (ones(36,1)*(eig_2_r+eig_3_r-2*eig_1_r)).*E1_x_E1 + ...
        (ones(36,1)*(eig_1_r)).*(E23_x_E1 + E1_x_E23) );
    % computation of the consistent tangent operators (36*n_int)
    Sder1_r=(ones(36,1)*(sigma_1_r-sigma_3_r)).*EIG_1_r;
    Sder2_r=VOL(:)*lame_r;
    Sder3_r=repmat(shear_r,36,1).*( 2*E1_x_E1 + E23_x_E23 );
    D_phi_r=repmat(shear_r,6,1).*(2*repmat(1+sin_phi_r,6,1).*Eig_1_r-repmat(1-sin_phi_r,6,1).*Eig_23_r)...
        +2*iota*(lame_r.*sin_phi_r);
    D_psi_r=repmat(shear_r,6,1).*(2*repmat(1+sin_psi_r,6,1).*Eig_1_r-repmat(1-sin_psi_r,6,1).*Eig_23_r)...
        +2*iota*(lame_r.*sin_psi_r);
    Sder4_r = ...
        [ D_psi_r(1,:).*D_phi_r(1,:); D_psi_r(2,:).*D_phi_r(1,:); D_psi_r(3,:).*D_phi_r(1,:); D_psi_r(4,:).*D_phi_r(1,:); D_psi_r(5,:).*D_phi_r(1,:); D_psi_r(6,:).*D_phi_r(1,:)
        D_psi_r(1,:).*D_phi_r(2,:); D_psi_r(2,:).*D_phi_r(2,:); D_psi_r(3,:).*D_phi_r(2,:); D_psi_r(4,:).*D_phi_r(2,:); D_psi_r(5,:).*D_phi_r(2,:); D_psi_r(6,:).*D_phi_r(2,:)
        D_psi_r(1,:).*D_phi_r(3,:); D_psi_r(2,:).*D_phi_r(3,:); D_psi_r(3,:).*D_phi_r(3,:); D_psi_r(4,:).*D_phi_r(3,:); D_psi_r(5,:).*D_phi_r(3,:); D_psi_r(6,:).*D_phi_r(3,:)
        D_psi_r(1,:).*D_phi_r(4,:); D_psi_r(2,:).*D_phi_r(4,:); D_psi_r(3,:).*D_phi_r(4,:); D_psi_r(4,:).*D_phi_r(4,:); D_psi_r(5,:).*D_phi_r(4,:); D_psi_r(6,:).*D_phi_r(4,:)
        D_psi_r(1,:).*D_phi_r(5,:); D_psi_r(2,:).*D_phi_r(5,:); D_psi_r(3,:).*D_phi_r(5,:); D_psi_r(4,:).*D_phi_r(5,:); D_psi_r(5,:).*D_phi_r(5,:); D_psi_r(6,:).*D_phi_r(5,:)
        D_psi_r(1,:).*D_phi_r(6,:); D_psi_r(2,:).*D_phi_r(6,:); D_psi_r(3,:).*D_phi_r(6,:); D_psi_r(4,:).*D_phi_r(6,:); D_psi_r(5,:).*D_phi_r(6,:); D_psi_r(6,:).*D_phi_r(6,:) ]./repmat(denom_r(test_r),36,1);
    DS(:,test_r)=Sder1_r+Sder2_r+Sder3_r-Sder4_r;

    % return to the apex
    DS(:,test_a)=zeros(36,nt_a);

end % (if nargout)

end % (function)
