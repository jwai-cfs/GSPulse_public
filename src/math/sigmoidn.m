function y = sigmoidn(x, x0, x1, y0, y1)
% =========================================================================
% Description: 
% Return a smooth transition (scaled sigmoid) that starts at (x0,y0) and
% ends at (x1,y1). Points outside of the x-range are held constant. 
%
% Example: 
% x0 = 0.2;
% x1 = 0.3;
% y0 = 0;
% y1 = 1;
% x = linspace(x0-0.1,x1+0.1)
% y = sigmoidn(x, x0, x1, y0, y1)
% plot(x,y,[x0 x1], [y0 y1], 'o', 'markerfacecolor', 'r', 'markersize', 12)
%
% =========================================================================

scale = 10 / (x1-x0);
shift = (x1+x0) / 2;
y = 1 ./ (1 + exp(-scale .* (x-shift)));
y(x<x0) = 0;
y(x>x1) = 1;
y = y0 + y*(y1-y0);