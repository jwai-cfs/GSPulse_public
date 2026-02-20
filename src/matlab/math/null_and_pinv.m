function [nullA, pinvA] = null_and_pinv(A)
% find nullspace and inv of A
% inspired by null.m and pinv.m, but more efficient since only one svd call

[m,n] = size(A);
[U,S,V] = svd(A);
s = diag(S);
tol = eps*10;
r = nnz(s > tol);

nullA = V(:,r+1:n);
pinvA = V(:,1:r) ./ s(1:r)' * U(:,1:r)';

assert(r >= m, 'GSpulse problem is infeasible, equality constraints cannot be satisfied');