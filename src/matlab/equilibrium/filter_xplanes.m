function [r, z, ikeep] = filter_xplanes(rq, zq, rmaxis, zmaxis, rx, zx)
% *************************************************************************
% Description: 
%
% Exclude x-points that lie outside the domain formed by the x-point
% half-planes. See [F.M. Joret, 2015, Fusion Engr. Design, section 2.3 
% and eqn 15]. 
% 
% Inputs: 
%  rq - R of query points
%  zq - Z of query points
%  rmaxis - R of magnetic axis
%  zmaxis - Z of magnetic axis
%  rx - R of x-points
%  zx - Z of x-points
% 
% Outputs: 
%  r - the valid R of the query points
%  z - the valid Z of the query points
%  ikeep - index of points that were kept
% *************************************************************************
out = false(size(rq));

for i = 1:length(rx)
  out = out | (rq-rx(i))*(rx(i)-rmaxis) + (zq-zx(i))*(zx(i)-zmaxis) > 0;
end

ikeep = ~out;
r = rq(ikeep);
z = zq(ikeep);