function [A, Ax, Ay, Axx, Ayy, Axy] = bicubicHermiteMat(xg, yg, x0, y0)
% =========================================================================
% Description:
% see bicubicHermite.m, returns matrices for 2D interpolation. This use
% case turns the 2D-interpolation operation into a matrix multiplication.
%
% Let [z, zx, zy, zxx, zyy, zxy] be the interpolated output of using
% bicubicHermite.m, i.e
%
% [z, zx, zy, zxx, zyy, zxy] = bicubicHermite(xg, yg, zg, x0, y0)
%
% The matrix outputs of this function are computed such that if we produce
% them from
%
% [A, Ax, Ay, Axx, Ayy, Axy] = bicubicHermiteMat(xg, yg, x0, y0)
%
% then the interpolated values can be produced as:
%
% z = A * zg(:);
% zx = Ax * zg(:);
% zy = Ay * zg(:);
% zxx = Axx * zg(:);
% zyy = Ayy * zg(:);
% zxy = Axy * zg(:);
% =========================================================================

% initialize
xg = xg(:);
yg = yg(:);
x0 = x0(:);
y0 = y0(:);
nx = length(xg);
ny = length(yg);
npts = numel(x0);
dx = xg(2) - xg(1);
dy = yg(2) - yg(1);
A   = sparse(npts, ny*nx);
Ax  = sparse(npts, ny*nx);
Ay  = sparse(npts, ny*nx);
Axx = sparse(npts, ny*nx);
Ayy = sparse(npts, ny*nx);
Axy = sparse(npts, ny*nx);

% grid point weighting matrix for cubic convolution (Keys, IEEE 1981)
mx = 1/2 * [0 2 0 0; -1 0 1 0; 2 -5 4 -1; -1 3 -3 1];

for k = 1:npts

  % if input point is nan then propagate the nan forward to the outputs
  if isnan(x0(k)) || isnan(y0(k)) || isinf(x0(k)) || isinf(y0(k))
    A(k,1) = nan;
    Ax(k,1) = nan;
    Ay(k,1) = nan;
    Axx(k,1) = nan;
    Ayy(k,1) = nan;
    Axy(k,1) = nan;
    continue;
  end

  % find the neigboring points (2 in each direction) around the query pt
  ix = find(xg < x0(k), 1, 'last');
  iy = find(yg < y0(k), 1, 'last');

  % check if inbounds
  inbounds = ~isempty(ix) && ~isempty(iy) && ...
    iy >= 2 && iy <= ny-2 && ix>=2 && ix <= nx-2;
  if ~inbounds
    A(k,1) = nan;
    Ax(k,1) = nan;
    Ay(k,1) = nan;
    Axx(k,1) = nan;
    Ayy(k,1) = nan;
    Axy(k,1) = nan;
    continue;
  end

  ineighbors = [ny*(ix-2) + [iy-1 iy iy+1 iy+2],
    ny*(ix-1) + [iy-1 iy iy+1 iy+2],
    ny*(ix+0) + [iy-1 iy iy+1 iy+2],
    ny*(ix+1) + [iy-1 iy iy+1 iy+2]];

  % normalize the x and y intervals
  tx = (x0(k) - xg(ix))/dx;
  ty = (y0(k) - yg(iy))/dy;

  % form the interpolation matrices
  by0 = [1 ty ty^2 ty^3] * mx;
  by1 = [0 1 2*ty 3*ty^2]/dy*mx;
  by2 = [0 0 2 6*ty]/dy^2*mx;

  bx0 = kron(eye(4), [1 tx tx^2 tx^3] * mx);
  bx1 = kron(eye(4), [0 1 2*tx 3*tx^2]/dx * mx);
  bx2 = kron(eye(4), [0 0 2 6*tx]/dx^2*mx);

  A(k,ineighbors)   = by0 * bx0;
  Ax(k,ineighbors)  = by0 * bx1;
  Ay(k,ineighbors)  = by1 * bx0;
  Axx(k,ineighbors) = by0 * bx2;
  Ayy(k,ineighbors) = by2 * bx0;
  Axy(k,ineighbors) = by1 * bx1;
end