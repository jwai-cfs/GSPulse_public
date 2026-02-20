x = 1:0.02:3;
y = -1.5: 0.02: 1.5;
[xgg,ygg] = meshgrid(x,y);
z = cos(5*xgg).*sin(ygg);
dx_smooth_span = 0.8;
dy_smooth_span = 0.3;
do_pad = true;
dx_pad = 0.8;
dy_pad = 0.5;
pad_method = 'nearest';
z_smooth = smooth2d(x, y, z, dx_smooth_span, dy_smooth_span, do_pad, ...
  dx_pad, dy_pad, pad_method);

% figure
% mesh(x,y,z)
% hold on
% mesh(x,y,z_smooth)
% shg
