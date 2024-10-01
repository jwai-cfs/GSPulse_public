function vq = interp1hold(x,v,xq,varargin)
% Description: 
% Interpolates and hold end point values
%
% (With interp1, the only way to hold endpoint values is to use the
% 'nearest' method, which also affects the interior.) 
%
% This has only been robustly tested for scalar signals. 
% 
% Input/Outputs: has the same call signature as matlab built-in interp1.m
%
% Example: 
% x = 1:1:10;
% v = sin(x);
% xq = -1:0.1:12;
% vq = interp1hold(x, v, xq, 'spline'); % spline interpolation for interior
% plot(xq,vq,x,v,'.','markersize', 40)

% edge case where there is only one sample point
if (length(x) == 1) && (x == xq)
  vq = v; 
else
  xq1 = xq(xq < min(x));
  xq2 = xq(xq >= min(x) & xq <= max(x));
  xq3 = xq(xq > max(x));
  
  if isempty(xq1)
    vq1 = [];
  else
    vq1 = interp1(x,v,xq1,'nearest','extrap');
  end
  vq2 = interp1(x,v,xq2, varargin{:});
  if isempty(xq3)
    vq3 = []; 
  else
    vq3 = interp1(x,v,xq3, 'nearest','extrap');
  end

  
  sz1 = size(vq1);
  sz2 = size(vq2);
  sz3 = size(vq3);
  
  dim = find((sz1 ~= sz2) & (sz2 ~= sz3));
  if isempty(dim), dim = 1; end
  vq = cat(dim(1), vq1, vq2, vq3);
end