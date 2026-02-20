% reference: https://en.wikipedia.org/wiki/Shoelace_formula
function A = polygonarea(x,y)
  
  x = [x(:); x(1)];
  y = [y(:); y(1)];
  n = length(x);

  i1 = 2:n;
  i2 = 1:(n-1);

  A = 0.5 * sum(x(i2).*y(i1) - x(i1).*y(i2));
end
