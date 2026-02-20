function make_gif_from_soln(soln, gsp_inputs, tok_data_struct, filename)
% =========================================================================
% Description: 
%   makes a gif of the equilibrium evolution from the GSPulse solution
%
% Inputs: 
%  soln - GSPulse solution struct
%  gsp_inputs - the GSPulse inputs read from config files, also see
%               run_pulse.m
%  filename - filename to save gif
% 
% Outputs: 
%  writes gif to file
%
% =========================================================================

times = soln.t;
N = length(times);

fig = figure;
fig.Position = [486 250 731 616];
for i = 1:N

  clf

  ax1 = axes();
  ax1.Position = [0.1300 0.1100 0.3347 0.8150];
  plot_shape(times(i), times, gsp_inputs.shapes, soln.eqs, gsp_inputs.tok)
  grid on
  ax1.XLim = [1.2 2.7];
  ax1.YLim = [-1.7 1.7];
  
  ax2 = copyobj(ax1, fig);
  ax2.Position = [0.5703 0.5838 0.3347 0.3412];
  ax2.XLim = [1.5 1.9];
  ax2.YLim =  [1.25 1.65];
  ax2.Title.Visible = 'off';
  grid on

  ax3 = copyobj(ax1, fig);
  ax3.Position = [0.5703 0.1100 0.3347 0.3412];
  ax3.XLim = [1.5 1.9];
  ax3.YLim = [-1.65 -1.25];
  ax3.Title.Visible = 'off';
  grid on
  % sgtitle(['                                             Pulse ' num2str(shot)], 'fontsize', 24)

  % write gif
  drawnow
  if i == 1
    gif(filename, 'overwrite', true, 'DelayTime', 0.05)
  else
    gif
  end
end

% plot shape targets
function plot_shape(t, times, shapes, eqs, tok)

  [~,i] = min(abs(t-times));
  ref = structts2struct(shapes, fieldnames(shapes), t);
  cla
  hold on  
  plot([nan nan], [nan nan], 'b', 'linewidth', 2)
  plot([nan nan], [nan nan], 'r', 'linewidth', 2)
  plot(ref.rx, ref.zx, 'xb', 'linewidth', 2, 'markersize', 14)
  try
    [rx, zx] = isoflux_xpFinder(tok.rg, tok.zg, eqs{i}.psizr, ref.rx, ref.zx);
    plot(rx, zx, 'xr', 'linewidth', 2, 'markersize', 14)  
  catch
  end  
  if isfield(ref, 'rstrike')
    scatter(ref.rstrike, ref.zstrike, 20, 'b', 'filled')
  end
  plot_eq(eqs{i}, tok, 'r', 'linewidth', 1.5)
  scatter(ref.rb, ref.zb, 20, 'b', 'filled')
  scatter(ref.rtouch, ref.ztouch, 100, 'db', 'filled')
  plot_lim(tok, 'k', 'linewidth', 2)
  plot_coils(tok_data_struct, [], 0.3)
  str = sprintf('Equilibrium %d: time=%.3f', i, t);
  title(str, 'fontsize', 14)
  legend('Target', 'Actual', 'fontsize', 14)
  xlabel('R [m]', 'fontsize', 12)
  ylabel('Z [m]', 'fontsize', 12)
  drawnow
end
end