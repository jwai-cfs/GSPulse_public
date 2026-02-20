function [rx, zx, psix, is_opoint] = isoflux_xpFinder(rg, zg, psizr, rx0, zx0)
%
% ISOFLUX_XPFINDER
%
%   Given the magnetic flux from a plasma equilibrium, determine the (r,z)
%   coordinates of a magnetic null given an initial guess (rx0,zx0) for the
%   null position.
%
% USAGE: isoflux_xpFinder.m
%
% INPUTS:
%
%   psizr...matrix with dimensions (nz x nr) containing the magnetic flux
%           at nz vertical by nr radial grid points
%
%   rx0.....initial guess for the null radial coordinate
%
%   zx0.....initial guess for the null vertical coordinate
%
%   rg......array containing the nr radial grid points
%
%   zg......array containing the nz vertical grid points
%
% OUTPUTS: 
%
%   rx......radial position of the null
%
%   zx......vertical position of the null
%
%   psix....magnetic flux at the null
%
% METHOD: This script uses Newton's method to locate the null (Br = Bz = 0)
%         and bicubic Hermite splines to interpolate on the grid.
%                       
%
%..........................................................................
e_relax = 0.3;
% Grid point weighting matrix for cubic convolution (Keys, IEEE 1981)

mx = 1/2 * [0 2 0 0; -1 0 1 0; 2 -5 4 -1; -1 3 -3 1];

% Convert rg and zg to column vectors if necessary

if size(rg,1) == 1
    rg = rg';
end
if size(zg,1) == 1
    zg = zg';
end

dr = rg(2)-rg(1);
dz = zg(2)-zg(1);

% Convert psizr to (nz x nr) if necessary

nr = length(rg);
nz = length(zg);

if size(psizr,1) == 1 || size(psizr,2) == 1
    psizr = reshape(psizr, nz, nr);
end

ii  = 100;  % maximum of 20 iterations to find the null
brzmax = inf;

rx = rx0;
zx = zx0;

while ii > 0 && brzmax > 1e-10
    
     ii = ii-1;
    
    % Locate four grid points in each direction around the query point
    
    ir = find(rg < rx, 1, 'last');
    iz = find(zg < zx, 1, 'last');
   
    ii1 = nz*(ir-2) + [iz-1 iz iz+1 iz+2];
    ii2 = nz*(ir-1) + [iz-1 iz iz+1 iz+2];
    ii3 = nz*(ir+0) + [iz-1 iz iz+1 iz+2];
    ii4 = nz*(ir+1) + [iz-1 iz iz+1 iz+2];
    
    F = [psizr(ii1(1)) psizr(ii1(2)) psizr(ii1(3)) psizr(ii1(4)); ...
         psizr(ii2(1)) psizr(ii2(2)) psizr(ii2(3)) psizr(ii2(4)); ...
         psizr(ii3(1)) psizr(ii3(2)) psizr(ii3(3)) psizr(ii3(4)); ...
         psizr(ii4(1)) psizr(ii4(2)) psizr(ii4(3)) psizr(ii4(4))  ...
    ];

    % Normalize the r and z intervals
     
    tr = (rx - rg(ir))/dr;
    tz = (zx - zg(iz))/dz;

    % Interpolate to find psi and its derivatives
     
    b0 = [1 tr tr^2 tr^3]*mx*F(:,1);
    b1 = [1 tr tr^2 tr^3]*mx*F(:,2);
    b2 = [1 tr tr^2 tr^3]*mx*F(:,3);
    b3 = [1 tr tr^2 tr^3]*mx*F(:,4);
    
    b0_r = ([0 1 2*tr 3*tr^2]/dr)*mx*F(:,1);
    b1_r = ([0 1 2*tr 3*tr^2]/dr)*mx*F(:,2);
    b2_r = ([0 1 2*tr 3*tr^2]/dr)*mx*F(:,3);
    b3_r = ([0 1 2*tr 3*tr^2]/dr)*mx*F(:,4);
    
    b0_rr = ([0 0 2 6*tr]/dr^2)*mx*F(:,1);
    b1_rr = ([0 0 2 6*tr]/dr^2)*mx*F(:,2);
    b2_rr = ([0 0 2 6*tr]/dr^2)*mx*F(:,3);
    b3_rr = ([0 0 2 6*tr]/dr^2)*mx*F(:,4);
     
    psix = [1 tz tz^2 tz^3]*mx*[b0 b1 b2 b3]';
    
    psi_r = [1 tz tz^2 tz^3]*mx*[b0_r b1_r b2_r b3_r]';
    psi_z = ([0 1 2*tz 3*tz^2]/dz)*mx*[b0 b1 b2 b3]';
    
    psi_rr = [1 tz tz^2 tz^3]*mx*[b0_rr b1_rr b2_rr b3_rr]';
    psi_zz = ([0 0 2 6*tz]/dz^2)*mx*[b0 b1 b2 b3]';
    psi_rz = ([0 1 2*tz 3*tz^2]/dz)*mx*[b0_r b1_r b2_r b3_r]';
        
    delta = -inv([psi_rr psi_rz; psi_rz psi_zz])*[psi_r; psi_z];
    rx = max(min(rx + e_relax * delta(1),rg(end-2)),rg(3));
    zx = max(min(zx + e_relax * delta(2),zg(end-2)),zg(3));
    
% scatter(rx,zx,'k','filled')
    brzmax = max(abs([psi_r,psi_z]));
    
end

is_opoint = 0;
% look at local curvature to determine if o or x-point
e = eig([psi_rr psi_rz; psi_rz psi_zz]);
if prod(sign(e)) ~= -1
  is_opoint = 1;
end

end
