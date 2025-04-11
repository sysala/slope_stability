% ************************************************************************

function   [S,DS]= constitutive_problem_2D(E,c_bar,sin_phi,shear,bulk,lame)
% =========================================================================
%
% The aim of this function is to construct constitutive and consistent 
% tangent operators at integration points 1,2,...,n_int. These operators
% are related to a simplified static version of the elastic-perfectly
% plastic model with the Mohr-Coulomb yield criterion.
%
% Input data:
%  E     - current strain tensor, size(E)=(3,n_int)
%  c_bar,sin_phi,shear,bulk,lame - material parameters at integration points
%
% Output data:
%  S      - projection of E at integration points, size(S)=(4,n_int)
%  DS     - derivative of the projection at integr. points,
%           size(DS)=(9,n_plast)
%  n_plast- number of the integration points with plastic response
%
% =========================================================================                                   

% 
% Linear constitutive operators
%
  IDENT = diag([1, 1, 1/2, 1]) ;         % identity operator  
  IOTA  = [1;1;0;1] ;                    % unit second tensor
  VOL   = IOTA*IOTA' ;                   % volumetric operator  
  Vol=VOL(1:3,1:3);
  Ident = diag([1, 1, 1/2]) ;
  Elast = Vol(:)*lame+2*Ident(:)*shear;  % reduced elastic operator

%
% Trial strain:
%   E_tr   - trial strain tensors, size(E_tr)=(4,n_int)
%
  n_int=size(E,2); % number of integration points
  E_tr=[E;zeros(1,n_int)];                     

%
% Inordered eigenvalues of the trial strain  
%
  I1=E_tr(1,:)+E_tr(2,:);
  I2=sqrt((E_tr(1,:)-E_tr(2,:)).^2+E_tr(3,:).^2);
  eig0_1=(I1+I2)/2 ;
  eig0_2=(I1-I2)/2 ;        
  eig0_3=E_tr(4,:) ;  

%  
% Inordered eigenprojections of the trial strain (1st der.)
%
  Eig0_1 = zeros(4,n_int);
  test1=(I2==0);
  Eig0_1(1:3,~test1) = [E_tr(1,~test1)-eig0_2(~test1)
                        E_tr(2,~test1)-eig0_2(~test1)
                        E_tr(3,~test1)/2]./...
                       (ones(3,1)*I2(~test1));  
  Eig0_1(1:2,test1) = ones(2,length(eig0_1(test1)));               
  Eig0_2 = [ones(2,n_int); zeros(2,n_int)]-Eig0_1;
  Eig0_3 = [zeros(3,n_int) ; ones(1,n_int)];

%  
% Inordered second derivatives of the eigenvalues
%
  EIG0_1 = zeros(9,n_int);
  EIG0_2 = zeros(9,n_int);
  EIG0_3 = zeros(9,n_int);   
  EIG0_1(:,~test1)=...
  [1-Eig0_1(1,~test1).*Eig0_1(1,~test1)-Eig0_2(1,~test1).*Eig0_2(1,~test1)
    -Eig0_1(2,~test1).*Eig0_1(1,~test1)-Eig0_2(2,~test1).*Eig0_2(1,~test1)
    -Eig0_1(3,~test1).*Eig0_1(1,~test1)-Eig0_2(3,~test1).*Eig0_2(1,~test1)
    -Eig0_1(1,~test1).*Eig0_1(2,~test1)-Eig0_2(1,~test1).*Eig0_2(2,~test1)
   1-Eig0_1(2,~test1).*Eig0_1(2,~test1)-Eig0_2(2,~test1).*Eig0_2(2,~test1)
    -Eig0_1(3,~test1).*Eig0_1(2,~test1)-Eig0_2(3,~test1).*Eig0_2(2,~test1)
    -Eig0_1(1,~test1).*Eig0_1(3,~test1)-Eig0_2(1,~test1).*Eig0_2(3,~test1)
    -Eig0_1(2,~test1).*Eig0_1(3,~test1)-Eig0_2(2,~test1).*Eig0_2(3,~test1)
  1/2-Eig0_1(3,~test1).*Eig0_1(3,~test1)-Eig0_2(3,~test1).*Eig0_2(3,~test1)...
  ]./(ones(9,1)*I2(~test1));
  EIG0_2(:,~test1)=-EIG0_1(:,~test1);

%  
% Reordering of eigenvalues and their derivatives  
%
  
  eig_1 = eig0_1; 
  eig_2 = eig0_2;   
  eig_3 = eig0_3;   % ordered eigenvalues of the trial strain
  Eig_1 = Eig0_1;
  Eig_2 = Eig0_2;
  Eig_3 = Eig0_3;   % eigenprojections of the trial strain (1st der.)
  EIG_1 = EIG0_1;
  EIG_2 = EIG0_2;
  EIG_3 = EIG0_3;   % the second derivatives of the eigenvalues
       
  test2=(eig0_1>=eig0_3)&(eig0_3>eig0_2);
  eig_2(test2)=eig0_3(test2);
  eig_3(test2)=eig0_2(test2);
  Eig_2(:,test2)=Eig0_3(:,test2);
  Eig_3(:,test2)=Eig0_2(:,test2);
  EIG_2(:,test2)=EIG0_3(:,test2);
  EIG_3(:,test2)=EIG0_2(:,test2);
  
  test3=(eig0_3>eig0_1);
  eig_1(test3)=eig0_3(test3);
  eig_2(test3)=eig0_1(test3);
  eig_3(test3)=eig0_2(test3);   % reordered eigenvalues of the trial strain
  Eig_1(:,test3)=Eig0_3(:,test3);
  Eig_2(:,test3)=Eig0_1(:,test3);
  Eig_3(:,test3)=Eig0_2(:,test3);
  EIG_1(:,test3)=EIG0_3(:,test3);
  EIG_2(:,test3)=EIG0_1(:,test3);
  EIG_3(:,test3)=EIG0_2(:,test3); % reordered derivatives
  
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
  S = zeros(4,n_int);
%   lambda = zeros(1,n_int);
%   sigma_1 = zeros(1,n_int); 
%   sigma_2 = zeros(1,n_int);   
%   sigma_3 = zeros(1,n_int);   % principal stresses
  
  % elastic response
  test_el=(f_tr<=0);
  lame_el=lame(test_el);
  shear_el=shear(test_el);
  nt_el=length(lame_el);
  S(:,test_el)=  repmat(lame_el,4,1).*(VOL*E_tr(:,test_el))+...
               2*repmat(shear_el,4,1).*(IDENT*E_tr(:,test_el)); 

  % return to the smooth portion of the yield surface
  test_s=(lambda_s<=min(gamma_sl,gamma_sr))&(~test_el);
  lambda_s=lambda_s(test_s);
  lame_s=lame(test_s);
  shear_s=shear(test_s);
  sin_phi_s=sin_phi(test_s);
  denom_s=denom_s(test_s);
  Eig_1_s=Eig_1(:,test_s);
  Eig_2_s=Eig_2(:,test_s);
  Eig_3_s=Eig_3(:,test_s);
  trace_E_s=trace_E(test_s);
  nt_s=length(lambda_s);
  sigma_1_s=lame_s.*trace_E_s+2*shear_s.*eig_1(test_s)-...
                  lambda_s.*(2*lame_s.*sin_phi_s+2*shear_s.*(1+sin_phi_s));
  sigma_2_s=lame_s.*trace_E_s+2*shear_s.*eig_2(test_s)-...
                  lambda_s.*(2*lame_s.*sin_phi_s);
  sigma_3_s=lame_s.*trace_E_s+2*shear_s.*eig_3(test_s)-...
                  lambda_s.*(2*lame_s.*sin_phi_s-2*shear_s.*(1-sin_phi_s));  
  S(:,test_s)=(ones(4,1)*sigma_1_s).*Eig_1_s+...
              (ones(4,1)*sigma_2_s).*Eig_2_s+...
              (ones(4,1)*sigma_3_s).*Eig_3_s;
              
  % return to the left edge of the yield surface             
  test_l=(gamma_sl<gamma_sr)&(lambda_l>=gamma_sl)&(lambda_l<=gamma_la)&...
         (~(test_el|test_s));
  lambda_l=lambda_l(test_l);  
  lame_l=lame(test_l);
  shear_l=shear(test_l);
  sin_phi_l=sin_phi(test_l);
  denom_l=denom_l(test_l);
  Eig_1_l=Eig_1(:,test_l);
  Eig_2_l=Eig_2(:,test_l);
  Eig_3_l=Eig_3(:,test_l);
  trace_E_l=trace_E(test_l);
  nt_l=length(lambda_l);
  sigma_1_l=lame_l.*trace_E_l+...
                  shear_l.*(eig_1(test_l)+eig_2(test_l))-...
                  lambda_l.*(2*lame_l.*sin_phi_l+shear_l.*(1+sin_phi_l));
  sigma_3_l=lame_l.*trace_E_l+2*shear_l.*eig_3(test_l)-...
                  lambda_l.*(2*lame_l.*sin_phi_l-2*shear_l.*(1-sin_phi_l));
  S(:,test_l)=(ones(4,1)*sigma_1_l).*(Eig_1_l+Eig_2_l)+...
              (ones(4,1)*sigma_3_l).*Eig_3_l;   
          
  % return to the right edge of the yield surface             
  test_r=(gamma_sl>gamma_sr)&(lambda_r>=gamma_sr)&(lambda_r<=gamma_ra)&...
         (~(test_el|test_s));
  lambda_r=lambda_r(test_r);  
  lame_r=lame(test_r);
  shear_r=shear(test_r);
  sin_phi_r=sin_phi(test_r);
  denom_r=denom_r(test_r);
  Eig_1_r=Eig_1(:,test_r);
  Eig_2_r=Eig_2(:,test_r);
  Eig_3_r=Eig_3(:,test_r);
  trace_E_r=trace_E(test_r);
  nt_r=length(lambda_r);
  sigma_1_r=lame_r.*trace_E_r+2*shear_r.*eig_1(test_r)-...
                  lambda_r.*(2*lame_r.*sin_phi_r+2*shear_r.*(1+sin_phi_r));
  sigma_3_r=lame_r.*trace_E_r+shear_r.*(eig_2(test_r)+eig_3(test_r))-...
                  lambda_r.*(2*lame_r.*sin_phi_r-shear_r.*(1-sin_phi_r));
  S(:,test_r)=(ones(4,1)*sigma_1_r).*Eig_1_r+...                                
              (ones(4,1)*sigma_3_r).*(Eig_2_r+Eig_3_r);        
              
  % return to the apex of the yield surface
  test_a=~(test_el|test_s|test_l|test_r);
  lambda_a=lambda_a(test_a);  
  nt_a=length(lambda_a);
  if (n_int~=nt_el+nt_s+nt_l+nt_r+nt_a), warning('number of elements!'), end
  sigma_1_a=c_bar(test_a)./(2*sin_phi(test_a));
  S(:,test_a)=IOTA*sigma_1_a;
     
%
% determination of the tangential stiffness array DS
%   
  if nargout>1
   % initialization
   DS = zeros(9,n_int);
   
   % elastic response
   DS(:,test_el)=Elast(:,test_el);
   
   % return to the smooth portion of the yield surface
   mat1_s=(ones(9,1)*sigma_1_s).*EIG_1(:,test_s)+...
          (ones(9,1)*sigma_2_s).*EIG_2(:,test_s)+...
          (ones(9,1)*sigma_3_s).*EIG_3(:,test_s);
   mat2_s=Vol(:)*lame_s;
   mat3_s=2*[shear_s.*Eig_1_s(1,:).*Eig_1_s(1,:)
             shear_s.*Eig_1_s(2,:).*Eig_1_s(1,:)
             shear_s.*Eig_1_s(3,:).*Eig_1_s(1,:)
             shear_s.*Eig_1_s(1,:).*Eig_1_s(2,:)
             shear_s.*Eig_1_s(2,:).*Eig_1_s(2,:)
             shear_s.*Eig_1_s(3,:).*Eig_1_s(2,:)
             shear_s.*Eig_1_s(1,:).*Eig_1_s(3,:)
             shear_s.*Eig_1_s(2,:).*Eig_1_s(3,:)
             shear_s.*Eig_1_s(3,:).*Eig_1_s(3,:) ];
   mat4_s=2*[shear_s.*Eig_2_s(1,:).*Eig_2_s(1,:)
             shear_s.*Eig_2_s(2,:).*Eig_2_s(1,:)
             shear_s.*Eig_2_s(3,:).*Eig_2_s(1,:)
             shear_s.*Eig_2_s(1,:).*Eig_2_s(2,:)
             shear_s.*Eig_2_s(2,:).*Eig_2_s(2,:)
             shear_s.*Eig_2_s(3,:).*Eig_2_s(2,:)
             shear_s.*Eig_2_s(1,:).*Eig_2_s(3,:)
             shear_s.*Eig_2_s(2,:).*Eig_2_s(3,:)
             shear_s.*Eig_2_s(3,:).*Eig_2_s(3,:) ];             
   mat5_s=2*[shear_s.*Eig_3_s(1,:).*Eig_3_s(1,:)
             shear_s.*Eig_3_s(2,:).*Eig_3_s(1,:)
             shear_s.*Eig_3_s(3,:).*Eig_3_s(1,:)
             shear_s.*Eig_3_s(1,:).*Eig_3_s(2,:)
             shear_s.*Eig_3_s(2,:).*Eig_3_s(2,:)
             shear_s.*Eig_3_s(3,:).*Eig_3_s(2,:)
             shear_s.*Eig_3_s(1,:).*Eig_3_s(3,:)
             shear_s.*Eig_3_s(2,:).*Eig_3_s(3,:)
             shear_s.*Eig_3_s(3,:).*Eig_3_s(3,:) ];           
   Eig_6=repmat(2*shear_s.*(1+sin_phi_s),3,1).*Eig_1_s(1:3,:)-...
         repmat(2*shear_s.*(1-sin_phi_s),3,1).*Eig_3_s(1:3,:)+...
         kron([1;1;0],2*lame_s.*sin_phi_s);
   mat6_s=      [ Eig_6(1,:).*Eig_6(1,:)
                  Eig_6(2,:).*Eig_6(1,:)
                  Eig_6(3,:).*Eig_6(1,:)
                  Eig_6(1,:).*Eig_6(2,:)
                  Eig_6(2,:).*Eig_6(2,:)
                  Eig_6(3,:).*Eig_6(2,:)
                  Eig_6(1,:).*Eig_6(3,:)
                  Eig_6(2,:).*Eig_6(3,:)
                  Eig_6(3,:).*Eig_6(3,:) ]./repmat(denom_s,9,1);
   DS(:,test_s)=mat1_s+mat2_s+mat3_s+mat4_s+mat5_s-mat6_s;
  
   % return to the left edge of the yield surface  
   Eig_12_l=Eig_1_l(1:3,:)+Eig_2_l(1:3,:);
   EIG_12=EIG_1(:,test_l)+EIG_2(:,test_l);
   mat1_l=(ones(9,1)*sigma_1_l).*EIG_12+...
          (ones(9,1)*sigma_3_l).*EIG_3(:,test_l);
   mat2_l=Vol(:)*lame_l;
   mat3_l= [shear_l.*Eig_12_l(1,:).*Eig_12_l(1,:)
            shear_l.*Eig_12_l(2,:).*Eig_12_l(1,:)
            shear_l.*Eig_12_l(3,:).*Eig_12_l(1,:)
            shear_l.*Eig_12_l(1,:).*Eig_12_l(2,:)
            shear_l.*Eig_12_l(2,:).*Eig_12_l(2,:)
            shear_l.*Eig_12_l(3,:).*Eig_12_l(2,:)
            shear_l.*Eig_12_l(1,:).*Eig_12_l(3,:)
            shear_l.*Eig_12_l(2,:).*Eig_12_l(3,:)
            shear_l.*Eig_12_l(3,:).*Eig_12_l(3,:)]; 
   mat5_l=2*[shear_l.*Eig_3_l(1,:).*Eig_3_l(1,:)
             shear_l.*Eig_3_l(2,:).*Eig_3_l(1,:)
             shear_l.*Eig_3_l(3,:).*Eig_3_l(1,:)
             shear_l.*Eig_3_l(1,:).*Eig_3_l(2,:)
             shear_l.*Eig_3_l(2,:).*Eig_3_l(2,:)
             shear_l.*Eig_3_l(3,:).*Eig_3_l(2,:)
             shear_l.*Eig_3_l(1,:).*Eig_3_l(3,:)
             shear_l.*Eig_3_l(2,:).*Eig_3_l(3,:)
             shear_l.*Eig_3_l(3,:).*Eig_3_l(3,:) ];           
   Eig_6= repmat(shear_l.*(1+sin_phi_l),3,1).*Eig_12_l-...
          repmat(2*shear_l.*(1-sin_phi_l),3,1).*Eig_3_l(1:3,:)+...
          kron([1;1;0],2*lame_l.*sin_phi_l);
   mat6_l=       [Eig_6(1,:).*Eig_6(1,:)
                  Eig_6(2,:).*Eig_6(1,:)
                  Eig_6(3,:).*Eig_6(1,:)
                  Eig_6(1,:).*Eig_6(2,:)
                  Eig_6(2,:).*Eig_6(2,:)
                  Eig_6(3,:).*Eig_6(2,:)
                  Eig_6(1,:).*Eig_6(3,:)
                  Eig_6(2,:).*Eig_6(3,:)
                  Eig_6(3,:).*Eig_6(3,:) ]./repmat(denom_l,9,1);
   DS(:,test_l)=mat1_l+mat2_l+mat3_l+mat5_l-mat6_l;
  
   % return to the right edge of the yield surface  
   Eig_23_r=Eig_2_r(1:3,:)+Eig_3_r(1:3,:);
   EIG_23=EIG_2(:,test_r)+EIG_3(:,test_r);
   mat1_r=(ones(9,1)*sigma_1_r).*EIG_1(:,test_r)+...
          (ones(9,1)*sigma_3_r).*EIG_23;
   mat2_r=Vol(:)*lame_r;
   mat3_r=2*[shear_r.*Eig_1_r(1,:).*Eig_1_r(1,:)
             shear_r.*Eig_1_r(2,:).*Eig_1_r(1,:)
             shear_r.*Eig_1_r(3,:).*Eig_1_r(1,:)
             shear_r.*Eig_1_r(1,:).*Eig_1_r(2,:)
             shear_r.*Eig_1_r(2,:).*Eig_1_r(2,:)
             shear_r.*Eig_1_r(3,:).*Eig_1_r(2,:)
             shear_r.*Eig_1_r(1,:).*Eig_1_r(3,:)
             shear_r.*Eig_1_r(2,:).*Eig_1_r(3,:)
             shear_r.*Eig_1_r(3,:).*Eig_1_r(3,:) ];     
   mat5_r= [shear_r.*Eig_23_r(1,:).*Eig_23_r(1,:)
            shear_r.*Eig_23_r(2,:).*Eig_23_r(1,:)
            shear_r.*Eig_23_r(3,:).*Eig_23_r(1,:)
            shear_r.*Eig_23_r(1,:).*Eig_23_r(2,:)
            shear_r.*Eig_23_r(2,:).*Eig_23_r(2,:)
            shear_r.*Eig_23_r(3,:).*Eig_23_r(2,:)
            shear_r.*Eig_23_r(1,:).*Eig_23_r(3,:)
            shear_r.*Eig_23_r(2,:).*Eig_23_r(3,:)
            shear_r.*Eig_23_r(3,:).*Eig_23_r(3,:) ];               
   Eig_6=repmat(2*shear_r.*(1+sin_phi_r),3,1).*Eig_1_r(1:3,:)-...
         repmat(shear_r.*(1-sin_phi_r),3,1).*Eig_23_r+...
         kron([1;1;0],2*lame_r.*sin_phi_r);
   mat6_r=       [Eig_6(1,:).*Eig_6(1,:)
                  Eig_6(2,:).*Eig_6(1,:)
                  Eig_6(3,:).*Eig_6(1,:)
                  Eig_6(1,:).*Eig_6(2,:)
                  Eig_6(2,:).*Eig_6(2,:)
                  Eig_6(3,:).*Eig_6(2,:)
                  Eig_6(1,:).*Eig_6(3,:)
                  Eig_6(2,:).*Eig_6(3,:)
                  Eig_6(3,:).*Eig_6(3,:) ]./repmat(denom_r,9,1);
   DS(:,test_r)=mat1_r+mat2_r+mat3_r+mat5_r-mat6_r;
  
   % return to the apex of the yield surface  
   DS(:,test_a)=zeros(9,nt_a);
  
  end % if nargout>1

 end
