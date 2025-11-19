function [pw, grad_p, mater_sat]=seepage_problem_3D...
           (coord,elem,Q_w,pw_D,grho,conduct0,HatP,DHatP1,DHatP2,DHatP3,WF)

  
  % number of nodes, elements, etc.
  n_e=size(elem,2);           % number of elements
  dim = size(coord, 1);         % Spatial dimension (2D or 3D).
  n_n=size(coord,2);          % number of nodes
  n_q=length(WF);             % number of quadratic points
  n_int = n_e*n_q ;           % total number of integrations points 
  n_p=size(elem,1);           % number of nodes per element

  % penalty parameter dependent on the element size
  eps=SEEPAGE.penalty_parameters_3D(coord,elem);

  % values of the penalization parameter at integration points
  eps_int=kron(eps,ones(1,n_q));

  % assembling of auxilliary matrices
  [Bw,C,weight]=SEEPAGE.auxiliary_matrices_3D...
                                 (elem,coord,HatP,DHatP1,DHatP2,DHatP3,WF);  
 
  % assembling of the constitutive matrix D, size(D)=(3*n_int, 3*n_int)
  wc=weight.*conduct0; 
  vD=[1;0;0;0;1;0;0;0;1]*wc; % size(vD)=(9, n_int) 
  AUX=reshape(1:3*n_int,3,n_int);
  iD=repmat(AUX,3,1); 
  jD=kron(AUX,ones(3,1));
  D=sparse(iD,jD,vD);
  
  % flow stiffness matrix: size(K_D)=(n_n, n_n)
  K_D = Bw'*D*Bw;     
  
  % the right-hand side vector
  % wc3=[wc; wc; wc];
  % q1=Bw*pw_D'; 
  % f=-Bw'*(wc3(:).*q1);

  wc2=repmat(wc,dim,1);
  q1=Bw*pw_D'; q2=Bw*coord(2,:)';
  q3=q1+grho*q2;
  f=-Bw'*(wc2(:).*q3);
  
  % initialization of the pressure field
  pw_0=zeros(1,n_n);
  pw_0(Q_w)=K_D(Q_w,Q_w)\f(Q_w);
  pw_init=pw_0+pw_D;

  % Newton's solver
  it_max=50; % maximal number of iterations
  tol=1e-10; % relative tolerance for the Newton solver
  pw=SEEPAGE.newton_flow(pw_init,conduct0,Q_w,weight,Bw,C,K_D,wc,...
            elem,coord,HatP,WF,eps_int,grho,it_max,tol);   

  % Remaining output arrays
  grad_p=reshape(Bw*pw',3,n_int);
  pw_e=reshape(pw(elem(:)),n_p,n_e);
  pw_int=sum(repmat(HatP,1,n_e).*kron(pw_e,ones(1,n_q)));
  if n_q>1
    int_pw_e=sum(reshape(pw_int.*weight,n_q,n_e));
    int_e=sum(reshape(weight,n_q,n_e));
    pw_aver_e=int_pw_e./int_e;
  else
    pw_aver_e=pw_int;
  end
  mater_sat=false(1,n_e);
  mater_sat(pw_aver_e>=0.1*eps)=1;    
                   
end