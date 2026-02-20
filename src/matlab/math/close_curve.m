function [x, y] = close_curve(x, y, tol)
% =========================================================================
% Description: 
%
% closes the curve defined by x and y, only if x and y dont already make a
% closed curve. x and y are considered an open curve if the distance
% between (x(1), y(1)) and (x(end), y(end)) is greater than tol. 
%
% Inputs: 
%  x - x coordinates of curve
%  y - y coordinates of curve
%  tol - numerical tolerance for checking if curve is alread closed
%
% Outputs: 
%  (x,y) - closed curve
%
% =========================================================================

if ~exist('tol', 'var') || isempty(tol)
  tol = sqrt(eps);
end

if sqrt((x(1)-x(end))^2 + (y(1)-y(end))^2) > tol
  x(end+1) = x(1);
  y(end+1) = y(1);
end