% function [xc, yc] = polygoncentroid(x,y)
% reference: https://mathworld.wolfram.com/PolygonCentroid.html
function [xc, yc] = polygoncentroid(x,y)

n = length(x);
xb = [x(:); x(1)];
yb = [y(:); y(1)];

i = 1:n;
j = 2:n+1;

A = polygonarea(x,y);
crossprod = (xb(i).*yb(j) - xb(j).*yb(i));
xc = sum((xb(i) + xb(j)) .* crossprod) / (6*A);
yc = sum((yb(i) + yb(j)) .* crossprod) / (6*A);