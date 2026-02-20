function [yq, BM] = spline_basis(x,y,xq)
% =========================================================================
% Description: 
%  Performs a 1-D cubic spline interpolation, and returns an equivalent
%  splining basis matrix. 
% 
%  yq = spline_basis(x, y, xq) returns the same output as
%  yq = interp1(x, y, x1, 'spline')
% 
%  The splining basis matrix BM is returned such that:
%  yq = BM*y
% 
% Inputs:
%   x - x grid for the input data
%   y - y values for the input data
%   xq - query points for the interpolator
%
% Outputs:
%   yq - interpolated values
%   BM - splining basis matrix
%
% Reference: https://en.wikipedia.org/wiki/Cubic_Hermite_spline
%
% Example:
% x = linspace(1, 5.3, 10)';
% y = sin(pi*x);
% xq = linspace(1, 5.3, 100)';
% [~, BM] = spline_basis(x,y,xq);
% yq = BM*y;
% 
% figure
% hold on
% scatter(x, y, 60, 'r', 'filled')
% plot(xq, sin(pi*xq), '--r')
% plot(xq, yq, 'b')
% legend('Interpolation points', 'True', 'Estimate', 'fontsize', 16)

% Restrictions:
%  - requires uniform spacing in x
%  - all of xq must lie within x interval (endpoints inclusive)
%
% =========================================================================
if numel(x) ~= numel(y), warning('x and y sizes do not match'); end

dx = mean(diff(x));
x0 = [1 1:length(x) length(x)]';  % integer basis for x
nx0 = length(x0);
y0 = [y(1); y(:); y(end)];
B = zeros(length(xq), length(y0));

for i = 1:length(xq)

  x1 = 1 + (xq(i) - x(1)) / dx;  
  n = floor(x1);
  u = x1 - n;
  
  % edge case, point lies at end of x-interval
  if abs(n+u+2 - nx0) < sqrt(eps)   
    b = 1;
    idx = length(y0);  

  % most cases    
  else                  
    b = 0.5 * [-u^3 + 2*u^2 - u; ...
      3*u^3 - 5*u^2 + 2; ...
      -3*u^3 + 4*u^2 + u; ...
      u^3 - u^2];

    idx = n + [-1 0 1 2] + 1; 
  end
  
  B(i,idx) = b;

end

BM = [B(:,1) + B(:,2),  B(:,3:end-2), B(:,end-1) + B(:,end)];
yq = BM*y;

end