function plot_coils(tok, c, alpha)
% =========================================================================
% Description: 
%  plot the coil geometries
%
% Inputs: 
%  tok - tokamak geometry object
%  c - color for which to plot coils, e.g. rgb vector or string like 'r'
%  alpha - transparency value for the coils
%
% Outputs: 
%  None
% =========================================================================

fcdata = tok.fcdata;
ecdata = [];
if isfield(tok, 'ecdata') && ~isempty(tok.ecdata)
  ecdata = tok.ecdata;
  xmin = min(ecdata(2,:) - ecdata(4,:)/2);
  xmax = max(ecdata(2,:) + ecdata(4,:)/2);
  ymin = min(ecdata(1,:) - ecdata(3,:)/2);
  ymax = max(ecdata(1,:) + ecdata(3,:)/2);
  xcenter = (xmin + xmax)/2;
  ycenter = (ymin + ymax)/2;
  dx = xmax - xmin;
  dy = ymax - ymin;
  ecdata = [ycenter xcenter dy dx 0 0]';  
end
ccdata = [ecdata fcdata];
nc = size(ccdata,2);


% housekeeping
if ~exist('c','var') || isempty(c)
  c = [1 0.75 0];
end

if ischar(c)
  c = bitget(find('krgybmcw'==c)-1,1:3);
end

if ~exist('alpha','var') || isempty(alpha)
  alpha = 1;
end

if size(c,1) ~= nc, c = ones(nc,1)*c; end
if isscalar(alpha), alpha = alpha*ones(nc,1); end


% plot coils
for i = 1:nc
  
%   y = ccdata(1,i) - ccdata(3,i)/2;
%   x = ccdata(2,i) - ccdata(4,i)/2;
%   dx = ccdata(4,i);
%   dy = ccdata(3,i);
%   position = [x y dx dy];
% 
%   h = rectangle('Position', position);
%   h.FaceColor = [c(i,:) alpha(i)];
%   h.EdgeColor = c(i,:) * min(1-alpha(i), 0.8);

  z = ccdata(1,i);
  r = ccdata(2,i);
  dz = ccdata(3,i);
  dr = ccdata(4,i);
  ac = ccdata(5,i);
  ac2 = ccdata(6,i);
  rgb = [0 0 0];
  
  if alpha(i) ~= 0
    h = plot_efit_region(z, r, dz, dr, ac, ac2, rgb);
    
    %   h = rectangle('Position', position);
    h.FaceColor = c(i,:);
    h.FaceAlpha = alpha(i);
    h.EdgeColor = c(i,:) * min(1-alpha(i), 0.8);
  end

end


% axis
xmax = max(ccdata(2,:));
xmin = min(ccdata(2,:));
ymax = max(ccdata(1,:));
ymin = min(ccdata(1,:));
axis equal
dx = (xmax-xmin)*0.1;
dy = (ymax-ymin)*0.1;

axis([xmin-dx xmax+dx ymin-dy ymax+dy])






































