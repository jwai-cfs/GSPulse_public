% normalizes each column of X to have vector norm 1
function [Xn, norms] = normc(X)
Xn = X;
n = size(X,2);
norms = nan(n,1);

for j = 1:n
  norms(j) = vecnorm(X(:,j));
  Xn(:,j) = X(:,j) / norms(j);
end
