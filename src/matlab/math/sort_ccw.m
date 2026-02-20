% sorts the points in (r,z) according to increasing counterclockwise angle
% that they make with (r0,z0). 

function [r,z] = sort_ccw(r,z,r0,z0)

if nargin == 2
  r0 = median(r);
  z0 = median(z);
end

[a,b] = dir_cosine(r0, z0, 0, r, z, r*0);
th = atan2(b,a);
[~,i] = sort(th);
r = r(i);
z = z(i);