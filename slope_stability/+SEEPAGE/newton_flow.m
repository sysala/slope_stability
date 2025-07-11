function pw=newton_flow(pw_init,conduct0,Q_w,weight,B,C,K_D,wc,...
                   elem,coord,HatP,WF,eps_int,grho,it_max,tol)

% Newton solver for unconfined seepage problem

% number of nodes, elements and integration points + print
  n_n=size(coord,2);            % number of nodes
  dim = size(coord, 1);         % Spatial dimension (2D or 3D).
  n_e=size(elem,2);             % number of elements
  n_q=length(WF);               % number of quadratic points
  n_int = n_e*n_q ;             % total number of integrations points 
  n_p=size(elem,1);             % number of nodes within one element
  
% initialization  
  pw=pw_init;
  it=0;
  
% iterations   
  while true
      
    it=it+1;
    % relative permeability and its derivative
    pw_e=reshape(pw(elem(:)),n_p,n_e);
    pw_int=sum(repmat(HatP,1,n_e).*kron(pw_e,ones(1,n_q)));
    perm_r=ones(1,n_int);
    perm_r_der=zeros(1,n_int);
    part1=(pw_int<eps_int)&(pw_int>0);
    part2=(pw_int<=0);
    perm_r(part1)=pw_int(part1)./eps_int(part1);
    perm_r(part2)=0;    
    perm_r_der(part1)=1./eps_int(part1);
       
    % assembling of the constitutive matrix E, size(E)=(dim*n_int, n_int)
    vE=repmat(conduct0.*perm_r_der.*weight,dim,1).*reshape(B*(grho*coord(2,:)'),dim,n_int); % size(vE)=(dim, n_int) 
    iE=reshape(1:dim*n_int,dim,n_int);
    jE=repmat(1:n_int,dim,1);
    E=sparse(iE,jE,vE); 
    
    % assembling of the stiffness matrix K
    K=K_D+B'*E*C;
    
    % assembling of the right-hand side vector
    kappa=repmat(perm_r,dim,1);
    wc2=repmat(wc,dim,1);
    q1=B*pw'; q2=B*coord(2,:)';
    q3=q1+grho*kappa(:).*q2;
    f=-B'*(wc2(:).*q3);
     
    % Newton's increment and next iteration
    dp=zeros(1,n_n);
    dp(Q_w)=K(Q_w,Q_w)\f(Q_w);
    pw=pw+dp;
    
    % stopping criteria
    crit=norm(dp)/norm(pw_init);
    if crit<tol
        fprintf(    'Newton solver converges: number of iteration=%d ',it); 
        fprintf(    'stopping criterion=%e ',crit);
        fprintf('\n'); 
        break
    end
    if it>it_max
        warning('Newton solver does not converge.')
        fprintf(    'stopping criterion=%e ',crit);
        fprintf('\n'); 
        break
    end
  end
    
end