function fig = summary_shape_plot(shapes, tok)
% =========================================================================
% Description: 
%   makes a slider plot of plasma target shapes
%
% Inputs: 
%  shapes - shapes struct, see help _define_shapes.m
%  tok - tokamak geometry struct, see help _define_tok.m
% 
% Outputs: 
%  fig - figure handle for the slider plot
% =========================================================================
tref = linspace(shapes.rb.Time(1), shapes.rb.Time(end), 200);

% create figure
fig = figure;
fig.Position = [754   322   425   611];
ax = axes(fig);
ax.Box = 'on';
ax.Position = [0.13 0.2 0.775 0.7];


% slider button
s = uicontrol(fig, 'style', 'slider');
s.Units = 'normalized';
s.Position = [0.15 0.05 0.6 0.05];
s.Min = 1;
s.Max = length(tref);
s.Value = 1;
s.SliderStep = [1 1] / (s.Max - s.Min);
s.Callback = {@sliderCallback, tref, shapes, tok};

plot_shape(s.Value, tref, shapes, tok)


% slider callback
function sliderCallback(src, event, tref, shapes, tok)
  t_idx = src.Value;
  plot_shape(t_idx, tref, shapes, tok)
end


% plot shape targets
function plot_shape(t_idx, tref, shapes, tok)

  t = tref(t_idx);

  % read parameters
  rb = interp1(shapes.rb.Time, shapes.rb.Data, t);
  zb = interp1(shapes.zb.Time, shapes.zb.Data, t);
  rx = interp1(shapes.rx.Time, shapes.rx.Data, t);
  zx = interp1(shapes.zx.Time, shapes.zx.Data, t);

  if isfield(shapes, 'rtouch')
    rtouch = interp1(shapes.rtouch.Time, shapes.rtouch.Data, t);
    ztouch = interp1(shapes.ztouch.Time, shapes.ztouch.Data, t);
  end

  if isfield(shapes, 'rstrike')
    rstrike = interp1(shapes.rstrike.Time, shapes.rstrike.Data, t);
    zstrike = interp1(shapes.zstrike.Time, shapes.zstrike.Data, t);
  end

  rb(end+1) = rb(1);
  zb(end+1) = zb(1);

  % plot parameters
  cla
  hold on
  title(['Time = ' num2str(t)], 'fontsize', 16)
  plot(tok.limdata(2,:), tok.limdata(1,:), 'k', 'linewidth', 1.5)
  % plot(rb, zb, 'b')
  scatter(rb, zb, 20, 'r', 'filled')
  plot(rx, zx, 'xb', 'linewidth', 4, 'markersize', 14)
  
  if isfield(shapes, 'rtouch')
    scatter(rtouch, ztouch, 100, 'db', 'filled')
  end

  if isfield(shapes, 'rstrike')
    rstrike = interp1(shapes.rstrike.Time, shapes.rstrike.Data, t);
    zstrike = interp1(shapes.zstrike.Time, shapes.zstrike.Data, t);
    scatter(rstrike, zstrike, 100, 'r', 'filled')
  end

  axis equal
  axis([min(tok.rg) max(tok.rg) min(tok.zg) max(tok.zg)])
  text(-0.25, -0.1, 'Drag slider to view shape targets.', ...
    'units', 'normalized', 'fontsize', 11)

  % axis([1.5509    1.8261   -1.6498   -1.1428])
  drawnow
end


end





































