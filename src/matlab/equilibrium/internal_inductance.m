function [Li, Wmag] = internal_inductance(ip, psizr, rb, zb, rg, zg)
% *************************************************************************
% Description: 
%  computes absolute (unnormalized) internal inductance
% 
% Inputs: 
%  ip - plasma current
%  psizr - flux on grid
%  (rb, zb) - R and Z of plasma boundary
%  rg - radial grid
%  zg - vertical grid
% 
% Outputs: 
%  Li - internal inductance
%  Wmag - stored magnetic energy
%
% *************************************************************************
[rgg,zgg] = meshgrid(rg,zg);
mu0 = pi*4e-7;
dr = mean(diff(rg));
dz = mean(diff(zg));
volumezr = 2*pi*rgg*dr*dz;

% [~,psi_r,psi_z] = bicubicHermite(rg, zg, psizr, rgg, zgg);

[psi_r, psi_z] = gradient(psizr, dr, dz);

Br = -1./(2*pi*rgg) .* psi_z;
Bz =  1./(2*pi*rgg) .* psi_r;

in = inpolygon(rgg, zgg, rb, zb);
Br(~in) = 0;
Bz(~in) = 0;

Wmag = 1 / (2*mu0) * sum(sum((Br.^2 + Bz.^2).*volumezr));  % magnetic field energy
Li = 2*Wmag / ip^2;










