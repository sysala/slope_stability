% ************************************************************************
%
%  constitutive_problem   -  solution of the constitutive 
%                            Drucker-Prager problem 
%
% ************************************************************************

function   Psi= potential_2D(E,c_bar,sin_phi,shear,bulk,lame)
% =========================================================================
%
% The aim of this function is to construct a potential to the constitutive 
% operator for the Mohr-Coulomb yield criterion and associative flow rule
% at integration points 1,2,...,n_int. 
%
% Input data:
%  E     - current strain tensor, size(E)=(3,n_int)
%  alpha - penalization parameter
%  shear,bulk,lame - parameters (constants)
%  c_bar,sin_phi - 1 x n_int arrays with heterogeneous parameters
%
% Output data:
%  Psi    - potential to E at integration points, size(S)=(1,n_int)
%
% =========================================================================                                   

%
% Trial strain:
%   E_tr   - trial strain tensors, size(E_tr)=(4,n_int)
%
  n_int=size(E,2); % number of integration points
  E_tr=[E;zeros(1,n_int)];                     

%
% Eigenvalues of the trial strain and their first and second derivatives
%
  
  I1=E_tr(1,:)+E_tr(2,:);
  I2=sqrt((E_tr(1,:)-E_tr(2,:)).^2+E_tr(3,:).^2);
  eig0_1=(I1+I2)/2 ;
  eig0_2=(I1-I2)/2 ;
  eig0_3=E_tr(4,:) ;  
 
% 
% Reordering of eigenvalues and their derivatives  
%
  
  eig_1 = eig0_1; 
  eig_2 = eig0_2;   
  eig_3 = eig0_3;   % ordered eigenvalues of the trial strain
       
  test2=(eig0_1>=eig0_3)&(eig0_3>eig0_2);
  eig_2(test2)=eig0_3(test2);
  eig_3(test2)=eig0_2(test2);
   
  test3=(eig0_3>eig0_1);
  eig_1(test3)=eig0_3(test3);
  eig_2(test3)=eig0_1(test3);
  eig_3(test3)=eig0_2(test3);   % reordered eigenvalues of the trial strain
   
%
% Critical values defining decision criteria
%
  trace_E=eig_1+eig_2+eig_3;        % trace of E_trial
  f_tr=2*shear.*((1+sin_phi).*eig_1-(1-sin_phi).*eig_3)+ ...
       2*lame.*sin_phi.*trace_E-c_bar; % test on admissibility
  gamma_sl=(eig_1-eig_2)./(1+sin_phi);
  gamma_sr=(eig_2-eig_3)./(1-sin_phi);
  gamma_la=(eig_1+eig_2-2*eig_3)./(3-sin_phi);
  gamma_ra=(2*eig_1-eig_2-eig_3)./(3+sin_phi);

%  
% Candidates on plastic multipliers
%
  denom_s=4*lame.*sin_phi.^2+2*shear.*(1+sin_phi).^2+2*shear.*(1-sin_phi).^2;
  denom_l=4*lame.*sin_phi.^2+  shear.*(1+sin_phi).^2+2*shear.*(1-sin_phi).^2;
  denom_r=4*lame.*sin_phi.^2+2*shear.*(1+sin_phi).^2+  shear.*(1-sin_phi).^2;
  denom_a=4*bulk.*sin_phi.^2;         % denominators for each type of return
  
  lambda_s=f_tr./denom_s ;
  lambda_l=(shear.*((1+sin_phi).*(eig_1+eig_2)-2*(1-sin_phi).*eig_3)+ ...
            2*lame.*sin_phi.*trace_E-c_bar)./denom_l ;
  lambda_r=(shear.*(2*(1+sin_phi).*eig_1-(1-sin_phi).*(eig_2+eig_3))+ ...
            2*lame.*sin_phi.*trace_E-c_bar)./denom_r ;
  lambda_a=(2*bulk.*sin_phi.*trace_E-c_bar)./denom_a ;

%
% Determination of the stress array S
%
  
  % initialization of unknowns variables
  Psi = zeros(1,n_int); 
  
  % elastic response
  test_el=(f_tr<=0);
  Psi(test_el)=0.5*lame(test_el).*trace_E(test_el).^2+shear(test_el).*(eig_1(test_el).^2+...
               eig_2(test_el).^2+eig_3(test_el).^2);
 
  % return to the smooth portion of the yield surface
  test_s=(lambda_s<=min(gamma_sl,gamma_sr))&(~test_el);
  Psi(test_s)=0.5*lame(test_s).*trace_E(test_s).^2+shear(test_s).*(eig_1(test_s).^2+...
              eig_2(test_s).^2+eig_3(test_s).^2)-...
              0.5*denom_s(test_s).*(lambda_s(test_s).^2);
                
  % return to the left edge of the yield surface             
  test_l=(gamma_sl<gamma_sr)&(lambda_l>=gamma_sl)&(lambda_l<=gamma_la)&...
         (~(test_el|test_s));
  Psi(test_l)=0.5*lame(test_l).*trace_E(test_l).^2+shear(test_l).*(eig_3(test_l).^2+...
              0.5*(eig_1(test_l)+eig_2(test_l)).^2)-...
              0.5*denom_l(test_l).*(lambda_l(test_l).^2);
            
  % return to the right edge of the yield surface             
  test_r=(gamma_sl>gamma_sr)&(lambda_r>=gamma_sr)&(lambda_r<=gamma_ra)&...
         (~(test_el|test_s));
  Psi(test_r)=0.5*lame(test_r).*trace_E(test_r).^2+shear(test_r).*(eig_1(test_r).^2+...
              0.5*(eig_2(test_r)+eig_3(test_r)).^2)-...
              0.5*denom_r(test_r).*(lambda_r(test_r).^2);  
              
  % return to the apex of the yield surface
  test_a=~(test_el|test_s|test_l|test_r);
  Psi(test_a)=0.5*bulk(test_a).*trace_E(test_a).^2-...
              0.5*denom_a(test_a).*(lambda_a(test_a).^2);
          
  end
