function [ri, zi] = seg_intersections(rb, zb, segs, plotit)
% =========================================================================
% Description: 
%  find the intersection points of the boundary with the control segments
%
% Inputs: 
%  (rb,zb) - plasma boundary
%  segs - control segments to find intersections with
%    The cols of segs are := [r_start r_end z_start z_end]
%    (r_start, z_start) should be the end point closest to the limiter, and
%    (r_end, z_end) should be the point closer to core. 
%
% Outputs:
%  (ri,zi) - intersect locations
% =========================================================================

if ~exist('plotit', 'var'), plotit = 0; end

% force to be loop
rb(end+1) = rb(1);
zb(end+1) = zb(1); 

for i = 1:size(segs,1)    
  [ri(i,1), zi(i,1)] = intersections(rb, zb, segs(i,1:2), segs(i,3:4));   
end


if plotit
  r0 = segs(:,1);
  rf = segs(:,2);
  z0 = segs(:,3);
  zf = segs(:,4);  
  plot(rb, zb)
  hold on
  plot([r0 rf]', [z0 zf]', 'color', [1 1 1] * 0.8, 'linewidth', 3)
  scatter(ri, zi, 60, 'b', 'filled')  
end


















