function Lz = calc_Lz(rg, zg, psizr, r0, z0, rBt, rmaxis, zmaxis, rlim, zlim, plotit)
% =========================================================================
% Description: 
%  post-processing step, compute the scrape-off-layer connection length 
%  of an equilibrium
% 
% Inputs:  
% rg - radial grid, 
% zg - vertical grid, 
% psizr - flux on grid, 
% (r0,z0) - seed location to start the trace, 
% (rmaxis,zmaxis) - magnetic axis
% (rlim,zlim) - the limiter contour
% plotit - flag for plotting
% 
% Outputs: 
% Lz - float, scrape off layer connection length
%
% Additional info: 
%   not really used for the GSPulse calculation, but needed by some
%   external that use GSPulse equilibria
% =========================================================================
opts = struct;
opts.step = 0;
opts.plotit = plotit;
opts.out_tol = 0.0001;
opts.check_every = 1;
opts.dt = 0.0002;

[rb, zb] = trace_contour(rg,zg,psizr,r0,z0,rmaxis,zmaxis,rlim,zlim,opts);
k = isnan(rb);
rb(k) = [];
zb(k) = [];

rc = (rb(1:end-1) + rb(2:end)) / 2;  % midpoint along boundary trace
zc = (zb(1:end-1) + zb(2:end)) / 2;

[~, psi_r, psi_z] = bicubicHermite(rg,zg,psizr,rc,zc);
Br =  1./(2*pi*rc).*psi_z;
Bz = -1./(2*pi*rc).*psi_r;
Bt = rBt ./ rc;

dr = diff(rb);
dz = diff(zb);
dt = dr .* Bt ./ Br;  % step length in toroidal direction, == dz .* Bt ./ Bz;
ds = sqrt(dr.^2 + dz.^2 + dt.^2); 
Lz = sum(ds);
