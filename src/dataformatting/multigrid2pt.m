function [ptvals, dptvalsdr, dptvalsdz] = multigrid2pt(rgrid, zgrid, gridvals, rpt, zpt)
% Description: 
% Interpolate multiple grids values to specific pts
%
% Inputs: rgrid - radial grid 
%         zgrid - vertical grid
%         gridvals - matrix of compressed gridded values, dimension [length(rgrid)*length(zgrid) x any]
%         rpt - r of points to be evaluated
%         zpt - z of points to be evaluated
% 
% Outputs: ptvals - matrix of pt values, dimension [npts x size(gridvals,2)]

nx = size(gridvals,2);
rpt = rpt(:);
zpt = zpt(:);
npt = length(rpt);
ptvals = zeros(npt,nx);
dptvalsdr = zeros(npt,nx);
dptvalsdz = zeros(npt,nx);


for i = 1:nx
  x = reshape(gridvals(:,i), length(zgrid), length(rgrid));
  [ptvals(:,i), dptvalsdr(:,i), dptvalsdz(:,i)] = bicubicHermite(...
    rgrid, zgrid, x, rpt, zpt);
end