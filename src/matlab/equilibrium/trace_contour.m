function [rb, zb] = trace_contour(rg,zg,psizr,r0,z0,rmaxis,zmaxis,rlim,zlim,opts)
% *************************************************************************
% Description: 
%   trace a contour of constant flux
% 
% Inputs:   
%   rg - radial grid, zg - vertical grid, psizr - flux on grid, (r0,z0) -
%   the seed point to start the tracing at, (rmaxis,zmaxis) - the magnetic
%   axis of the plasma, not required to actually be the magnetic axis but
%   it must be a point on the interior of the plasma, (rlim, zlim) -
%   description of the limiter, opts - if opts.plotit = true, makes a plot
%   of the equilibrium
%  
% Outputs:    
%  (rb,zb) - the traced boundary
%
% *************************************************************************

if ~exist('opts','var'),     opts = struct; end
if ~isfield(opts, 'plotit'), opts.plotit = 1; end
if ~isfield(opts, 'step'), opts.step = 0.005; end
if ~isfield(opts, 'dt'), opts.dt = 0.001; end   % parameter that affects step size
if ~isfield(opts, 'out_tol'), opts.out_tol = 0.01; end
if ~isfield(opts, 'check_every'), opts.check_every = 30; end
if ~isfield(opts, 'num_grad_desc_steps'), opts.num_grad_desc_steps = 2; end


% move initial condition slightly toward magnetic axis, improves robustness
vec = [rmaxis-r0; zmaxis-z0];
vec = vec/norm(vec);
r0 = r0 + opts.step*vec(1);
z0 = z0 + opts.step*vec(2);

% initialize
theta = 0;
N = 1e4;
rb = nan(N,1);
zb = nan(N,1);
rb(1) = r0;
zb(1) = z0;
psi0 = bicubicHermite(rg,zg,psizr,r0,z0);

for i = 1:N

  % orthogonal (gradient descent step) down to the psi=psi0 contour
  for j = 1:opts.num_grad_desc_steps
    [psi, psi_r, psi_z, psi_rr, psi_zz, psi_rz] = bicubicHermite(rg,zg,psizr,rb(i),zb(i));
    eps = 1;
    Jinv = pinv([psi_r psi_z]) * eps;
    rb(i) = rb(i) - Jinv(1) * (psi-psi0);
    zb(i) = zb(i) - Jinv(2) * (psi-psi0);
  end
  
  % parallel step along the psi0=constant contour
  f = [psi_z; -psi_r];
  gradf_f = [psi_rz*psi_z - psi_zz*psi_r; -psi_rr*psi_z + psi_rz*psi_r];
  drz = f*opts.dt + gradf_f*opts.dt^2/2;

  rb(i+1) = rb(i) + drz(1);
  zb(i+1) = zb(i) + drz(2);

  % measure change in angle of rotation, after 2pi rotation terminate
  u = [rb(i+1) 0 zb(i+1)] - [rmaxis 0 zmaxis];
  v = [rb(i) 0 zb(i)] - [rmaxis 0 zmaxis];
  theta = theta + atan2(norm(cross(u,v)),dot(u,v));
  tol_angle = 0;  
  if abs(theta) > 2*pi + tol_angle
    rb = rb(1:i);
    zb = zb(1:i);    
    break;
  end

  % if outside limiter terminate
  if mod(i,opts.check_every) == 0
    if ~inpolygon(rb(i), zb(i), rlim, zlim)
      d = sqrt(min((rlim-rb(i)).^2 + (zlim-zb(i)).^2));
      if d > opts.out_tol, break; end
    end
  end
end

if opts.plotit
  figure
  hold on
  contour(rg, zg, psizr, 20, 'color', [1 1 1]*0.8)
  plot(rlim, zlim, 'k')
  scatter(rb,zb,10,'b','filled')
  plot(rb,zb,'r','linewidth',2)
  axis equal
end
end