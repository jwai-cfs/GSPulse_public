function fig = shape_slider_plot(times, shapes, eqs, tok, L)
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
if isscalar(times)
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
  shape_fields = fieldnames(shapes);
  ref = structts2struct(shapes, shape_fields, t, 0);
  
  cla
  hold on  
  plot([nan nan], [nan nan], 'b', 'linewidth', 2)
  plot([nan nan], [nan nan], 'r', 'linewidth', 2)
  axis equal
  axis([min(tok.rg) max(tok.rg) min(tok.zg) max(tok.zg)])

  % plot target x-points
  for i = 1:4
    plot(ref.(['rx' num2str(i)]), ref.(['zx' num2str(i)]), 'xb', 'linewidth', 2, 'markersize', 14)
  end

  % plot target strike points
  for i = 1:8
    scatter(ref.(['rstrike' num2str(i)]), ref.(['zstrike' num2str(i)]), 20, 'b', 'filled')
  end

  % plot target strike points
  for i = 1:3
    scatter(ref.(['r_control_pt_ref' num2str(i)]), ref.(['z_control_pt_ref' num2str(i)]), 100, 'bs')
  end
  
  % plot contour lines along any found x-points
  [rA,zA,FA,dr2FA,dz2FA,drzFA,rX,zX,FX,dr2FX,dz2FX,drzFX,...
  rB,zB,FB,lB,lX,Opy,F0,F1,status,msg,id,dF0dFx,dF1dFx,ixI] = meqpdom(eqs{t_idx}.psizr,eqs{t_idx}.Ip,1,L);
  plot(rX, zX, 'xk', 'linewidth', 2, 'markersize', 14)

  psi_max = max(eqs{t_idx}.psizr(~tok.outside_vessel_mask(:)));
  psi_min = min(eqs{t_idx}.psizr(~tok.outside_vessel_mask(:)));

  % r = linspace(min(tok.rl), max(tok.rl), 8);
  % z = zeros(size(r));
  % levels = bicubicHermite(tok.rg, tok.zg, eqs{t_idx}.psizr, r, z);
  levels = linspace(psi_min, psi_max, 8);
  levels = [levels(:); FX(:)];
  contour(L.G.rx, L.G.zx, eqs{t_idx}.psizr, levels, 'color', [1 1 1] * 0.7, 'linewidth', 0.5)

  % plot equilibrium
  plot(tok.rl, tok.zl, 'k')
  if ~isempty(eqs{t_idx}.FB)
    contour(tok.rg, tok.zg, eqs{t_idx}.psizr, [1 1]*eqs{t_idx}.FB, 'r', 'linewidth', 1.5);  
  end

  % plot control points
  scatter(ref.cp_r, ref.cp_z, 20, 'b', 'filled')
  
  text(-0.25, -0.1, 'Drag slider to view equilibria', 'units', 'normalized', 'fontsize', 11)  
  text(1.05, -0.1, 'Enter time:', 'units', 'normalized', 'fontsize', 11)  
  str = sprintf('Equilibrium %d: time=%.3f', t_idx, t);
  title(str, 'fontsize', 14)
  legend('Target', 'Actual', 'fontsize', 14)
  grid on
  drawnow


end
end