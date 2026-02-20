function varargout = multigrid2pt(rgrid, zgrid, gridvals, rpt, zpt)
%
% [ptvals, dptvalsdr, dptvalsdz] = multigrid2pt(rgrid, zgrid, gridvals, rpt, zpt)
%
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
%          dptvalsdr,dptvalsdz: at desired points r and z derivatives


rpt = rpt(:);
zpt = zpt(:);

if size(gridvals,1) ~= length(zgrid) * length(rgrid) && ...
    size(gridvals,2) == length(zgrid) * length(rgrid)
  gridvals = gridvals';
end

x = reshape(gridvals, length(zgrid), length(rgrid),[]);
if nargout == 1
  [ptvals] = bicubicHermite(...
    rgrid, zgrid, x, rpt, zpt);
else
  [ptvals, dptvalsdr, dptvalsdz] = bicubicHermite(...
    rgrid, zgrid, x, rpt, zpt);
end

varargout{1} = ptvals;
if nargout > 1
 varargout{2} = dptvalsdr;
 varargout{3} = dptvalsdz;
end
end
