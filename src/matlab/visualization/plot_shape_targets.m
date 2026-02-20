function fig = plot_shape_targets(shapes, rl, zl)
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
tref = linspace(shapes.cp_r.Time(1), shapes.cp_r.Time(end), 200);

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
s.Callback = {@sliderCallback, tref, shapes, rl, zl};

plot_shape(s.Value, tref, shapes, rl, zl)


% slider callback
function sliderCallback(src, event, tref, shapes, rl, zl)
  t_idx = src.Value;
  plot_shape(t_idx, tref, shapes, rl, zl)
end


% plot shape targets
function plot_shape(t_idx, tref, shapes, rl, zl)

  t = tref(t_idx);
  
  cla
  hold on
  title(['Time = ' num2str(t)], 'fontsize', 16)
  plot(rl, zl, 'k', 'linewidth', 1.5)

  % read parameters
  cp_r = interp1(shapes.cp_r.Time, shapes.cp_r.Data, t);
  cp_z = interp1(shapes.cp_z.Time, shapes.cp_z.Data, t);
  scatter(cp_r, cp_z, 20, 'b', 'filled')

  for j = 1:4
    rx = interp1(shapes.(['rx' num2str(j)]).Time, shapes.(['rx' num2str(j)]).Data, t);
    zx = interp1(shapes.(['zx' num2str(j)]).Time, shapes.(['zx' num2str(j)]).Data, t);
    plot(rx, zx, 'xb', 'linewidth', 4, 'markersize', 14)  
  end
  
  for j = 1:8
    rstrike = interp1(shapes.(['rstrike' num2str(j)]).Time, shapes.(['rstrike' num2str(j)]).Data, t);
    zstrike = interp1(shapes.(['zstrike' num2str(j)]).Time, shapes.(['zstrike' num2str(j)]).Data, t);
    scatter(rstrike, zstrike, 50, 'b', 'filled')
  end
  
  for j = 1:3
    r = interp1(shapes.(['r_control_pt_ref' num2str(j)]).Time, shapes.(['r_control_pt_ref' num2str(j)]).Data, t);
    z = interp1(shapes.(['z_control_pt_ref' num2str(j)]).Time, shapes.(['z_control_pt_ref' num2str(j)]).Data, t);
    scatter(r, z, 150, 'bs', 'linewidth', 2)
  end
    
  axis equal
  
  text(-0.25, -0.1, 'Drag slider to view shape targets.', ...
    'units', 'normalized', 'fontsize', 11)

  % axis([1.5509    1.8261   -1.6498   -1.1428])
  drawnow
end


end





































