function [Br, Bz, Bp] = psizr2b(psizr, tok_data_struct, opts)
% *************************************************************************
% Description: 
%   extract magnetic fields on grid from flux distribution
% 
% Inputs: 
%   psizr - flux on grid
%   tok_data_struct - toksys geometry object
%
% Outputs: 
%   Br, Bz, Bp - radial, vertical and poloidal field magnitudes
%
% *************************************************************************
if ~exist('opts', 'var') || isempty(opts), opts = struct; end
if ~isfield(opts, 'plotit'), opts.plotit = 0; end

rg = tok_data_struct.rg;
zg = tok_data_struct.zg;
dr = mean(diff(rg));
dz = mean(diff(zg));
[rgg, zgg] = meshgrid(rg,zg);

% [~, psi_r, psi_z] = bicubicHermite(rg, zg, psizr, rgg, zgg);
[psi_r, psi_z] = gradient(psizr, dr, dz);

Br = -1./(2*pi*rgg)  .* psi_z;
Bz =  1 ./ (2*pi*rgg) .* psi_r;
Bp = sqrt(Br.^2 + Bz.^2);

if opts.plotit  
  contour(rg, zg, Br, 20)
end



