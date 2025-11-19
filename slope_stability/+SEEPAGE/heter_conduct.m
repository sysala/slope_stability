function conduct=heter_conduct(mater,n_q,k)
%--------------------------------------------------------------------------
% This function creates 1 x n_int array specifying the conductivity at
% each integration point.
%
%
% INPUT ARGUMENTS:
%   mater   - Array of material identifiers for each element. It is
%             assumed that these identifiers start from 0.
%   n_q     - Number of quadrature points per element.
%   k       - Array containing conductivities for each subdomain
%
% OUTPUT:
%   conduct - Array of k values for each integration point, size: (1, n_int)
%
%--------------------------------------------------------------------------
                                 
% numbers of elements and subdomain
  n_e=length(mater);

% Conductivity for each element
conduct = zeros(1, n_e);
for i = 1:n_e
    id = mater(i) + 1;
    conduct(i) = k(id);
end

% Conductivity for each integration point  
  conduct=kron(conduct,ones(1,n_q));  
              
end