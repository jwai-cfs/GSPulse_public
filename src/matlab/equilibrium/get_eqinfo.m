function [eq, warmstart_out] = get_eqinfo(rg, zg, psizr, rl, zl, opts)
% *************************************************************************
% Description: 
%  extract a lot of equilibrium info from the flux distribution
% 
% Inputs: 
%  rg - radial grid
%  zg - vertical grid
%  psizr - flux on the grid
%  rl - R of limiter
%  zl - Z of limiter
%  opts - options struct
% 
% Outputs: 
%  eq - struct with lots of equilibrium info such as the x-points, touch
%  points, and boundary
%
% *************************************************************************
% set default vars
if ~exist('opts', 'var'), opts = struct; end
default = struct;
default.plotit = 0;
default.nlim = 500; 
default.return_bry = 0;
default.warmstart = struct;
opts = copyfields(default, opts, fieldnames(opts), 1);

% grid info
dr = mean(diff(rg));
dz = mean(diff(zg));
nz = length(zg);
nr = length(rg);
psizr = reshape(psizr, nz, nr);
[rgg, zgg] = meshgrid(rg, zg);

% close and wrap the limiter, helps with finding extrema
[rl, zl] = close_curve(rl,zl);
rl = rl(:);
zl = zl(:);
rl = [rl; rl(2:4)];  
zl = [zl; zl(2:4)];

% compute warmstart data if not available
if ~isfield(opts.warmstart, 'rlimfine')
  [rlimfine, zlimfine] = interparc(rl, zl, opts.nlim, 1, 1);
  Ainterp = bicubicHermiteMat(rg, zg, rlimfine, zlimfine);  
  inside_lim_mask = inpolygon(rgg, zgg, rl, zl);
  opts.warmstart = variables2struct(rlimfine, zlimfine, Ainterp, inside_lim_mask);
end
warmstart_out = opts.warmstart;

% initialize search for nulls by selecting a very inclusive set of local
% extrema grid points. 
imin =  islocalmin(psizr, 'FlatSelection', 'all');
imax =  islocalmax(psizr, 'FlatSelection', 'all');
jmin =  islocalmin(psizr, 2, 'FlatSelection', 'all');
jmax =  islocalmax(psizr, 2, 'FlatSelection', 'all');

imin = get_neighbor_mask(imin, nr, nz);
imax = get_neighbor_mask(imax, nr, nz);
jmin = get_neighbor_mask(jmin, nr, nz);
jmax = get_neighbor_mask(jmax, nr, nz);

io = ((imin & jmin) | (imax & jmax)) & opts.warmstart.inside_lim_mask;  % index of possible o-points
ix = ((imin & jmax) | (imax & jmin)) & opts.warmstart.inside_lim_mask;  % index of possible x-points
ro = rgg(io);
zo = zgg(io);
rx = rgg(ix);
zx = zgg(ix);

% find unique nulls to within tolerance (pre-zoom)
tol = sqrt(dr^2+dz^2);
[~,i] = uniquetol([ro zo], tol, 'byrows', 1);
ro = ro(i);
zo = zo(i);
[~,i] = uniquetol([rx zx], tol, 'byrows', 1);
rx = rx(i);
zx = zx(i);

% zoom in to higher-than-grid resolution on x-points
psix = rx*nan;
is_opt_check = rx*nan;
for i = 1:length(rx)
  [rx(i), zx(i), psix(i), is_opt_check(i)] = isoflux_xpFinder(rg, zg, psizr, rx(i), zx(i));
end
is_opt_check = logical(is_opt_check);
rx(is_opt_check) = [];
zx(is_opt_check) = [];
psix(is_opt_check) = [];

% zoom in to higher-than-grid resolution on o-points
psio = ro*nan;
is_opt_check = ro*nan;
for i = 1:length(ro)
  [ro(i), zo(i), psio(i), is_opt_check(i)] = isoflux_xpFinder(rg, zg, psizr, ro(i), zo(i));
end

% Filter o-points to make sure still an o-point (didn't migrate to another
% null during the zoom). Helps eliminate spurious points from consideration 
% as magnetic axis. 
ro(~is_opt_check) = [];
zo(~is_opt_check) = [];
psio(~is_opt_check) = [];

% zooming in can move o-points outside limiter, exclude these from
% consideration as mag axis. 
[in, on] = inpolygon(ro, zo, rl, zl);
out = ~(in & ~on);
ro(out) = [];
zo(out) = [];
psio(out) = [];

% find unique nulls to within tolerance (post-zoom)
tol = 0.1 * sqrt(dr^2+dz^2);
[~,i] = uniquetol([ro zo], tol, 'byrows', 1);
ro = ro(i);
zo = zo(i);
psio = psio(i);

if length(rx) > 1
  [~,i] = uniquetol([rx zx], tol, 'byrows', 1);
  rx = rx(i);
  zx = zx(i);
  psix = psix(i);
end

% select the magnetic axis
if length(ro) < 1
  warning('No magnetic axis found')
elseif length(ro) == 1
  rmaxis = ro;
  zmaxis = zo;
  psimag = psio;
else
  warning('Multiple axes found. Selecting axis closest to center.');
  rc = (min(rl) + max(rl))/2;
  zc = (min(zl) + max(zl))/2;
  [~,i] = min((ro-rc).^2 + (zo-zc).^2);
  rmaxis = ro(i);
  zmaxis = zo(i);
  psimag = psio(i);
end

% determine curvature of magnetic axis (maxis_curvature=1 <==> local min)
[~, ~, ~, psi_rr, psi_zz] = bicubicHermite(rg, zg, psizr, rmaxis, zmaxis);
maxis_curvature = sign(sign(psi_rr) + sign(psi_zz));

% search for touch points by checking for local extrema along limiter flux
rlim_pts = opts.warmstart.rlimfine;
zlim_pts = opts.warmstart.zlimfine;
psilim_pts = opts.warmstart.Ainterp * psizr(:);
  
i = islocalmin(maxis_curvature*psilim_pts, 'FlatSelection', 'all');
rtouch = rlim_pts(i);
ztouch = zlim_pts(i);
psitouch = psilim_pts(i);

% Exclude xpts that lie outside the domain formed by the xpt semiplanes
[rx, zx, i] = filter_xplanes(rx, zx, rmaxis, zmaxis, rx, zx);
psix = psix(i);

% Exclude touchpoints that lie outside the domain from xpt semiplanes
[rtouch, ztouch, i] = filter_xplanes(rtouch, ztouch, rmaxis, zmaxis, rx, zx);
psitouch = psitouch(i);

% choose the remaining global extremum as the touch point
[~,i] = min(maxis_curvature*psitouch);
rtouch = rtouch(i);
ztouch = ztouch(i);
psitouch = psitouch(i);

% determine the boundary point from among the x-points and touch point 
rcand = [rx; rtouch];
zcand = [zx; ztouch];
psicand = [psix; psitouch];

[~,i] = min(maxis_curvature * psicand);
i = i(end);
rbdef = rcand(i);
zbdef = zcand(i);
psibry = psicand(i);
islimited = psibry == psitouch; 

% trace the boundary
if opts.return_bry
  trace_opts = struct;
  trace_opts.plotit = 0;
  trace_opts.dt = 0.001;
  [rbbbs, zbbbs] = trace_contour(rg,zg,psizr,rbdef,zbdef,rmaxis,zmaxis,rl,zl,trace_opts);
else
  rbbbs = [];
  zbbbs = [];
end

% define the plasma domain mask
psinzr = (psizr - psimag) / (psibry - psimag);
out = psinzr > 1 | psinzr < 0 | ~opts.warmstart.inside_lim_mask;
[~,~,inplasma] = filter_xplanes(rgg, zgg, rmaxis, zmaxis, rx, zx);
outside_plasma_domain_mask = out | ~inplasma;

% save data
eq = variables2struct(psizr, rbdef, zbdef, psibry, rmaxis, zmaxis, ...
  psimag, rx, zx, psix, rtouch, ztouch, psitouch, islimited, rbbbs, ...
  zbbbs, rg, zg, rl, zl, outside_plasma_domain_mask);

% plot results
if opts.plotit   
  clf
  hold on  
  contour(rg, zg, psizr, 50, 'color', [1 1 1]*0.8)
  scatter(ro, zo, 100, 'b', 'filled')
  scatter(rx, zx, 100, 'r', 'filled')
  scatter(rtouch, ztouch, 100, 'g', 'filled')
  scatter(rbdef, zbdef, 200, 'b', 'linewidth', 1.5)
  contour(rg, zg, psizr, [1 1]*psibry, 'r')
  plot(rl, zl, 'k')
  axis equal  
  drawnow
  shg
end