function out_mask = get_neighbor_mask(in_mask, nr, nz)
% *************************************************************************
% Description: 
%   Points are true in the out_mask if they are true in the in_mask, or if
%   they are a neighbor of a point that is true in the in_mask. A neighbor
%   is defined as the 8 points in a box immediately surrounding a point, or
%   the 4 points at (r+2,z), (r-2,z), (r,z+2), (r,z-2)
% 
% Inputs: 
%  in_mask - input mask
%  nr - number of radial grid points
%  nz - number of vertical grid points
% 
% Outputs: 
%  out_mask - output mask
% *************************************************************************
neighbors = [-2*nz 2*nz 2 -2 nz-1 nz nz+1 -1 1 -nz+1 -nz -nz-1];

k = find(in_mask) + neighbors;
k(k<1) = [];
k(k>nz*nr) = [];

out_mask = in_mask;
out_mask(k) = true;