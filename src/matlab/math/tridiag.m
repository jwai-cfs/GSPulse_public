function M = tridiag(v1, v2, v3, N)
% =========================================================================
% Description: 
% returns a tridiagonal matrix of size NxN with values v1 on the lower
% diagonal, v2 on the diagonal, and v3 on the upper diagonal. v1, v2, v3
% can be scalars or vectors of the appropriate size. 
% =========================================================================
if isscalar(v1), v1 = v1*ones(N-1,1); end
if isscalar(v2), v2 = v2*ones(N,1); end
if isscalar(v3), v3 = v3*ones(N-1,1); end

M = diag(v1,-1) + diag(v2) + diag(v3,1);