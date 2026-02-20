function z_smooth = smooth1d(x, z, dx_smooth_span, do_pad, dx_pad, pad_method)

% check inputs
if ~exist('do_pad', 'var'), do_pad = true; end
if ~exist('dx_pad', 'var'), dx_pad = dx_smooth_span; end
if ~exist('pad_method', 'var'), pad_method = 'nearest'; end

assert(isuniform(x), 'x grid must be uniform')

if dx_pad < dx_smooth_span
  warning(['Recommended to set pad size larger than smooth_span size, ' ...
    'to avoid edge distortion.'])
end

% interpret grid info
nx = length(x);
dx = mean(diff(x));

% span must odd, add 1 if even
nx_span = ceil(dx_smooth_span / dx);

if mod(nx_span,2)==0, nx_span = nx_span + 1; end

if ~do_pad
  z_smooth = apply_smoothing_kernel(z, nx_span);
else

  n_dx_pad = ceil(dx_pad / dx);
  x_pad = x(1) + dx * (-n_dx_pad:(nx+n_dx_pad-1));

  % Pad the array
  z_pad = interp1(x, z, x_pad, pad_method, 'extrap');

  % smooth
  z_pad_smooth = apply_smoothing_kernel(z_pad, nx_span);

  ix = n_dx_pad + (1:nx);
  z_smooth = z_pad_smooth(ix);
end
end

function z_smooth = apply_smoothing_kernel(z, nx_span)

  % initialize
  nx = length(z);
  z_smooth = nan(nx,1);
  
  % loop over grid
  for j = 1:nx
    idx_x = j + (-(nx_span-1)/2:(nx_span-1)/2); % x kernel indices
    idx_x(idx_x < 1) = [];  % truncate kernel to fit within grid min
    idx_x(idx_x > nx) = []; % truncate kernel to fit within grid max
  
    % take the mean
    z_smooth(j) = mean(z(idx_x));
  end
end
