function s = shape_analysis(r,z)
% =========================================================================
% Description: 
%  given a closed curve (r,z) representing the plasma boundary, measure 
%  various shaping parameters such as elongation, triangularity,
%  squareness, etc.
% 
% Inputs: 
%   (r,z) - input curve
%
% Outputs: 
%  s - struct with values set for (rsurf,zsurf) the geometric center,
%  aminor the minor radius, bminor the vertical minor radius, elong the
%  plasma elongation, aspectratio the aspect ratio, epsilon (1/aspect
%  ratio), triu/tril/tri the upper, lower, and total triangularity,
%  squo/squi/sqlo/sqli the squareness for the upper/lower and inner/outer
%  corners. 
%  
% =========================================================================

[r,z] = sort_ccw(r,z); 

% find inner, outer, upper, lower points
[ri,ii] = min(r);
[ro,io] = max(r);
[zu,iu] = max(z);
[zl,il] = min(z);
zi = z(ii);
zo = z(io);
ru = r(iu);
rl = r(il);

rsurf = (ro+ri)/2;
zsurf = (zu+zl)/2;
aminor = (ro-ri)/2;
bminor = (zu-zl)/2;
elong = bminor / aminor;
epsilon = aminor/rsurf;
aspectratio = 1/epsilon;
triu = (rsurf - ru) / aminor;
tril = (rsurf - rl) / aminor;
tri = (triu + tril) / 2;

% order matters for the squareness inputs 
% (outer/inner point should precede upper/lower point) 
squo = squareness(ro, zo, ru, zu, r, z);
sqlo = squareness(ro, zo, rl, zl, r, z);
squi = squareness(ri, zi, ru, zu, r, z);
sqli = squareness(ri, zi, rl, zl, r, z);

s = variables2struct(rsurf, zsurf, aminor, bminor, elong, triu, ...
  tril, tri, squo, sqlo, squi, sqli, ro, zo, ru, zu, ri, zi, rl, zl);

end

% squareness definition from: 
% https://iopscience.iop.org/article/10.1088/0741-3335/55/9/095009/meta
function sq = squareness(r1, z1, r2, z2, r, z)

  A = r1 - r2;
  B = z2 - z1;
  
  rellipse = linspace(r2, r1);
  zellipse = z1 + sign(B)*sqrt(B^2 - ((B/A)*(rellipse - r2)).^2);

  [rc, zc] = intersections([r2 r1], [z1 z2], rellipse, zellipse);
  [rd, zd] = intersections([r2 r1], [z1 z2], r, z);

  LOD = norm([rd - r2, zd - z1]);
  LOC = norm([rc - r2, zc - z1]);
  LCE = norm([rc - r1, zc - z2]);

  sq = (LOD - LOC) / LCE;

  % debug purposes
  if 0
    hold on
    plot(r,z)
    scatter([r1 r2], [z1 z2], 'filled')
    plot(rellipse, zellipse, '--k')
    scatter([rc rd], [zc zd], 20, 'dk', 'filled')
    axis equal
    axis([0 1.9 -1.2 1.2])
    title(sq)
  end

end