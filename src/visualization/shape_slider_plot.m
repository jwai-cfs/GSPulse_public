function fig = shape_slider_plot(times, shapes, eqs, tok)
% =========================================================================
% Description: 
%   makes a slider plot of the plasma equilibrium vs target shapes
%
% Inputs: 
%  times - vector of times corresponding to each of the equilibria
%  shapes - shapes struct, see help _define_shapes.m
%  eqs - cell array of equilibria
%  tok - tokamak geometry struct, see help _define_tok.m
% 
% Outputs: 
%  fig - figure handle for the slider plot
% =========================================================================

% create figure
fig = figure;
fig.Position = [754   322   425   611];
ax = axes(fig);
ax.Box = 'on';
ax.Position = [0.13 0.2 0.775 0.7];

% edge case where there is only 1 time
if length(times) == 1
  times = times * [1 1];
end

% slider button
s = uicontrol(fig, 'style', 'slider');
s.Units = 'normalized';
s.Position = [0.15 0.05 0.6 0.05];
s.Min = 1;
s.Max = length(times);
s.Value = 1;
s.SliderStep = [1 1] / (s.Max - s.Min);
s.Callback = {@sliderCallback, times, shapes, eqs, tok};

plot_shape(s.Value, times, shapes, eqs, tok)

% slider callback
function sliderCallback(src, event, times, shapes, eqs, tok)
  t_idx = src.Value;
  t_idx = round(t_idx);
  plot_shape(t_idx, times, shapes, eqs, tok)
end


% plot shape targets
function plot_shape(t_idx, times, shapes, eqs, tok)

  t = times(t_idx);
  ref = structts2struct(shapes, fieldnames(shapes), t);
  cla
  hold on  
  plot([nan nan], [nan nan], 'b', 'linewidth', 2)
  plot([nan nan], [nan nan], 'r', 'linewidth', 2)
  plot(ref.rx, ref.zx, 'xb', 'linewidth', 2, 'markersize', 14)
  try
    for j = 1:length(ref.rx)
      [rx, zx] = isoflux_xpFinder(tok.rg, tok.zg, eqs{t_idx}.psizr, ref.rx(j), ref.zx(j));
      plot(rx, zx, 'xr', 'linewidth', 2, 'markersize', 14)  
    end
  catch
  end  
  if isfield(ref, 'rstrike')
    scatter(ref.rstrike, ref.zstrike, 20, 'b', 'filled')
  end
  try
    plot_eq(eqs{t_idx}, tok, 'r', 'linewidth', 1.5)
  catch
  end
  plot_lim(tok)
  scatter(ref.rb, ref.zb, 20, 'b', 'filled')
  scatter(ref.rtouch, ref.ztouch, 100, 'db', 'filled')
  text(-0.25, -0.1, 'Drag slider to view equilibria', 'units', 'normalized', 'fontsize', 11)  
  text(1.05, -0.1, 'Enter time:', 'units', 'normalized', 'fontsize', 11)  
  str = sprintf('Equilibrium %d: time=%.3f', t_idx, t);
  title(str, 'fontsize', 14)
  legend('Target', 'Actual', 'fontsize', 14)
  grid on
  drawnow


end
end