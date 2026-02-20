function [rs, zs] = seg_strike_finder(rg, zg, psizr, psivals, segs, opts)
% *************************************************************************
% Description: 
% Use Newtons method to find 'strike points', locations along a segment
% that have a flux value equal to psivals. 
% 
% Inputs: 
%  rg - radial grid
%  zg - vertical grid
%  psizr - flux on grid
%  psivals - values of psi for which to find intersect locations
%  segs - segments, has dimensions [# segs x 4] where each row corresponds
%  to [r_start, r_end, z_start, z_end] of the segment
% 
% Outputs: 
%  (rs, zs) - the intersect locations for each of the psivals and segs
%
% *************************************************************************

% define default options
default.tol = 1e-8;
default.nmax = 100;
default.ds_tol = 0.001;
default.ds_max = 0.03;
default.plotit = 0;
if ~exist('opts', 'var'), opts = struct; end
opts = copyfields(default, opts, {}, 1);
if size(segs,2) ~= 4, segs = segs'; end

% initialize
nvals = length(psivals);
nsegs = size(segs,1);
rs = nan(nsegs, nvals);
zs = nan(nsegs, nvals);

% loop over segments and psivals to find intersection locations
for i = 1:nsegs

  % read seg geometry
  r0 = (segs(i,1) + segs(i,2)) / 2;
  z0 = (segs(i,3) + segs(i,4)) / 2;
  dr = segs(i,2) - segs(i,1);
  dz = segs(i,4) - segs(i,3);
   
  rmin = min(r0-dr/2, r0+dr/2);
  rmax = max(r0-dr/2, r0+dr/2);
  zmin = min(z0-dz/2, z0+dz/2);
  zmax = max(z0-dz/2, z0+dz/2);

  l = sqrt(dr^2 + dz^2);
  dr = dr / l;
  dz = dz / l; 

  
  for j = 1:nvals

    % initialize
    psival = psivals(j);
    psi = inf;
    r = r0;
    z = z0;

    % newtons method root finding along segment
    for k = 1:opts.nmax

      [psi, psi_r, psi_z] = bicubicHermite(rg, zg, psizr, r, z);
      dpsids = psi_r*dr + psi_z*dz;
      ds = (psival - psi) / dpsids;
      ds = sign(ds) * min(abs(ds), opts.ds_max);
      r = r + dr*ds;
      z = z + dz*ds;         

      % if meets tolerance for being an intersection, save and exit
      if abs(psival-psi) < opts.tol
        rs(i,j) = r;
        zs(i,j) = z;
        break
      end           
    end    
    % intersection is beyond the length of the segment, specify as nan
    if r > rmax+opts.ds_tol || r < rmin-opts.ds_tol || z > zmax+opts.ds_tol || z < zmin-opts.ds_tol
      rs(i,j) = nan;
      zs(i,j) = nan;  
    end
  end
end


% plottting
if opts.plotit
  if length(psivals) == 1, psivals = [1 1] * psivals; end
  hold on
  contour(rg, zg, psizr, 50, 'color', [1 1 1] * 0.95)
  contour(rg, zg, psizr, psivals, 'k')
  plot(segs(:,1:2)', segs(:,3:4)', 'b')
  scatter(rs, zs, 40, 'r', 'filled')
  axis([min(rg) max(rg) min(zg) max(zg)])
  axis equal
end



