% EXAMPLE: 
% x = 1:0.02:3;
% y = -1.5: 0.02: 1.5;
% [xgg,ygg] = meshgrid(x,y);
% z = cos(5*xgg).*sin(ygg);
% dx_smooth_span = 0.8;
% dy_smooth_span = 0.3;
% do_pad = true;
% dx_pad = 0.8;
% dy_pad = 0.5;
% pad_method = 'nearest';
% z_smooth = smooth2d(x, y, z, dx_smooth_span, dy_smooth_span, do_pad, ...
%   dx_pad, dy_pad, pad_method);
% figure
% mesh(x,y,z)
% hold on
% mesh(x,y,z_smooth)
% shg


function z_smooth = smooth2d(x, y, z, dx_smooth_span, dy_smooth_span, ...
  do_pad, dx_pad, dy_pad, pad_method)

% check inputs
if ~exist('do_pad', 'var'), do_pad = true; end
if ~exist('dx_pad', 'var'), dx_pad = dx_smooth_span; end
if ~exist('dy_pad', 'var'), dy_pad = dy_smooth_span; end
if ~exist('pad_method', 'var'), pad_method = 'nearest'; end

% assert(isuniform(x), 'x grid must be uniform')
% assert(isuniform(y), 'y grid must be uniform')

if dx_pad < dx_smooth_span || dy_pad < dy_smooth_span
  warning(['Recommended to set pad size larger than smooth_span size, ' ...
    'to avoid edge distortion.'])
end

% interpret grid info
nx = length(x);
ny = length(y);

dx = mean(diff(x));
dy = mean(diff(y));

nx_span = ceil(dx_smooth_span / dx);
ny_span = ceil(dy_smooth_span / dy);

if ~do_pad
  z_smooth = apply_smoothing_kernel(z, nx_span, ny_span);
else

  n_dx_pad = ceil(dx_pad / dx);
  n_dy_pad = ceil(dy_pad / dy);

  x_pad = x(1) + dx * (-n_dx_pad:(nx+n_dx_pad-1));
  y_pad = y(1) + dy * (-n_dy_pad:(ny+n_dy_pad-1));

  % Pad in 2D by using 2 subsequent 1D extrapolations. Unfortunately, no
  % simple built-in methods for 2D extrapolation. Its beneficial to pad
  % before smoothing so that the smoothing operation will not distort the 
  % edges. 
  z_pad = interp1(x, z', x_pad, pad_method, 'extrap');
  z_pad = interp1(y, z_pad',  y_pad, pad_method, 'extrap');

  % smooth
  z_pad_smooth = apply_smoothing_kernel(z_pad, nx_span, ny_span);
  
  ix = n_dx_pad + (1:nx);
  iy = n_dy_pad + (1:ny);
  z_smooth = z_pad_smooth(iy,ix);
end
end

function z_smooth = apply_smoothing_kernel(z, nx_span, ny_span)
  
  % span must odd, add 1 if even
  if mod(nx_span,2)==0, nx_span = nx_span + 1; end
  if mod(ny_span,2)==0, ny_span = ny_span + 1; end
  
  % initialize
  [ny,nx] = size(z);
  z_smooth = nan(ny,nx);
  
  % loop over grid
  for i = 1:ny
    idx_y = i + (-(ny_span-1)/2:(ny_span-1)/2);  % y kernel indices
    idx_y(idx_y < 1) = [];   % truncate kernel to fit within grid min
    idx_y(idx_y > ny) = [];  % truncate kernel to fit within grid max
  
    for j = 1:nx
      idx_x = j + (-(nx_span-1)/2:(nx_span-1)/2); % x kernel indices
      idx_x(idx_x < 1) = [];  % truncate kernel to fit within grid min
      idx_x(idx_x > nx) = []; % truncate kernel to fit within grid max
  
      % take the mean
      z_smooth(i,j) = mean(z(idx_y,idx_x), 'all');
    end
  end
end
