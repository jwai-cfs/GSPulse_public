function [rb, zb, rcp, zcp, rstrike, zstrike, rx, zx] = ...
  load_shape_file(shapefn, plotit)
% =========================================================================
% Description: 
%  loads a shape file (saved by the python shape GUI to a .json file) and
%  reads the shaping parameters into matlab
%
% Inputs: 
%  shapefn - the filename to load from
%  plotit - whether to plot the results
%  
% Outputs: 
%  (rb,zb) - the plasma boundary
%  (rcp, zcp) - the control points (intersection points of the boundary and
%  the control segments)
%  (rstrike, zstrike) - strike point locations
%  (rx, zx) - x point locations
%
% =========================================================================

if ~exist('plotit', 'var'), plotit = 0; end
d = load_json_dict(shapefn);

% full boundary
rb = d.shape_params.rb;
zb = d.shape_params.zb;
i = ~isnan(rb);
rb = rb(i);
zb = zb(i);

% boundary control points
rcp = d.shape_params.rcp;
zcp = d.shape_params.zcp;
i = ~isnan(rcp);
rcp = rcp(i);
zcp = zcp(i);

% strike points
for i = 1:8
  rfd = ['r' num2str(i)];
  zfd = ['z' num2str(i)];
  rstrike(i) = d.shape_params.(rfd);
  zstrike(i) = d.shape_params.(zfd);
end
i = ~isnan(rstrike);
rstrike = rstrike(i);
zstrike = zstrike(i);

% x-points
for i = 1:4
  rfd = ['rx' num2str(i)];
  zfd = ['zx' num2str(i)];
  if isfield(d.shape_params, rfd)
    rx(i) = d.shape_params.(rfd);
    zx(i) = d.shape_params.(zfd);
  end
end
i = ~isnan(rx);
rx = rx(i);
zx = zx(i);

% plot
if plotit     
  hold on
  plot(rl, zl, 'k', 'linewidth', 1.5)
  plot(rb, zb, 'r')
  scatter(rcp, zcp, 20, 'b')  
  axis equal
  grid on  
  scatter(rstrike, zstrike, 20, 'b')
  plot(rx, zx, 'rx', 'markersize', 12, 'linewidth', 3)
end