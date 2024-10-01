function [r,z] = shape_edit(r,z,s)
% =========================================================================
% Description:
%  given an input boundary curve (r,z), and the shaping parameters like
%  elongation/triangularity/squareness stored in s, edit the input curve to
%  have the properties of s.
%
% Inputs:
%   (r,z) - input curve
%   s - struct with shaping parameters
%
% Outputs:
%   (r,z) - the edited shape curve
%
% reference for squareness:
% https://iopscience.iop.org/article/10.1088/0741-3335/55/9/095009/meta
%
% The parameters in s that will edit the shape are: {rsurf, zsurf, aminor, 
% elong, triu, tril, squo, sqlo, squi, sqli}.  All other parameters are
% ignored during shape_edit.
% =========================================================================
[r,z] = sort_ccw(r,z);

s0 = shape_analysis(r,z);

% shape edits from (rsurf, zsurf)
r = r + s.rsurf - s0.rsurf;
z = z + s.zsurf - s0.zsurf;

% shape edits from aminor
r = s.rsurf + (r-s.rsurf) * s.aminor / s0.aminor;

% shape edits from elong
b0 = s0.aminor*s0.elong;
bminor = s.aminor * s.elong;
z = s.zsurf + (z-s.zsurf) * bminor / b0;

% shape edits from (triu, tril)
s0 = shape_analysis(r,z);

ru = s.rsurf - s.aminor * s.triu;
dru = ru - s0.ru;         % how much ru needs to move to match triu
iu = find(z > s0.zsurf);
f = interp1([s0.ri s0.ru s0.ro], [0 1 0], r(iu));
r(iu) = r(iu) + f * dru;  % movement is propto dru and distance from ri,ro

rl = s.rsurf - s.aminor * s.tril;
drl = rl - s0.rl;
il = find(z < s0.zsurf);
f = interp1([s0.ri s0.rl s0.ro], [0 1 0], r(il));
r(il) = r(il) + f * drl;

% shape edits from squareness
s0 = shape_analysis(r,z);

% order matters for the edit_squareness inputs
% (outer/inner point should precede upper/lower point)
[r1,z1] = edit_squareness(s0.ro, s0.zo, s0.ru, s0.zu, s0.squo, s.squo, r, z);
[r2,z2] = edit_squareness(s0.ri, s0.zi, s0.ru, s0.zu, s0.squi, s.squi, r, z);
[r3,z3] = edit_squareness(s0.ro, s0.zo, s0.rl, s0.zl, s0.sqlo, s.sqlo, r, z);
[r4,z4] = edit_squareness(s0.ri, s0.zi, s0.rl, s0.zl, s0.sqli, s.sqli, r, z);

r = [r1; r2; r3; r4];
z = [z1; z2; z3; z4];

[r,z] = sort_ccw(r,z);
r(end+1) = r(1);
z(end+1) = z(1);

end

function [r,z] = edit_squareness(r1,z1,r2,z2,sqinput,sqtarget,r,z)

bminor = z2 - z1;
aminor = r1 - r2;

% (x,y) is the (r,z) normalized to the quadrant 1 unit circle
x = (r-r2)./aminor;
y = (z-z1)./bminor;
i = x>=0 & y>=0;    % use only quadrant 1
x = x(i);
y = y(i);
[x,y] = sort_ccw(x,y);
th = cart2pol(x,y);

% curveA: the normalized input curve
curveA = variables2struct(x,y,th);

% curveB: the superellipse that matches input squareness, see ref
n = -log(2) / log(1/sqrt(2) + sqinput*(1-1/sqrt(2)));
x = linspace(1,0);
y =  (1 - x.^n) .^ (1/n);
th = cart2pol(x,y);
curveB = variables2struct(x,y,th);

% curveC: the superellipse that matches target squareness, see ref
n = -log(2) / log(1/sqrt(2) + sqtarget*(1-1/sqrt(2)));
x = linspace(1,0);
y =  (1 - x.^n) .^ (1/n);
th = cart2pol(x,y);
curveC = variables2struct(x,y,th);

% interpolate all curves according to angle
curveB.x = interp1(curveB.th, curveB.x, curveA.th);
curveB.y = interp1(curveB.th, curveB.y, curveA.th);
curveC.x = interp1(curveC.th, curveC.x, curveA.th);
curveC.y = interp1(curveC.th, curveC.y, curveA.th);

% shift input curveA by the amount that the superellipse shifted
dx = curveC.x - curveB.x;
dy = curveC.y - curveB.y;

% curveD: the normalized output curve
curveD.x = curveA.x + dx;
curveD.y = curveA.y + dy;

% denormalize
r = curveD.x*aminor + r2;
z = curveD.y*bminor + z1;
end