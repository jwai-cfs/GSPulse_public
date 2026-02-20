function fig = slider_plot_plasma_startup(soln, gsp_inputs)

% create figure
fig = figure;
fig.Position =  [918 388 903 531];
ax = axes(fig);
ax.Box = 'on';
ax.Position = [0.13 0.2 0.775 0.7];

% slider button
s = uicontrol(fig, 'style', 'slider');
s.Units = 'normalized';
s.Position = [0.15 0.05 0.6 0.05];
s.Min = 1;
s.Max = length(soln.t);
s.Value = 1;
s.SliderStep = [1 1] / (s.Max - s.Min);
s.Callback = {@sliderCallback, soln, gsp_inputs};

update_plots(s.Value, soln, gsp_inputs)

% slider callback
function sliderCallback(src, ~, soln, gsp_inputs)
  t_idx = src.Value;
  t_idx = round(t_idx);
  update_plots(t_idx, soln, gsp_inputs)
end

end

% update plots
function update_plots(t_idx, soln, gsp_inputs)
  
  subplot(131)
  plot_equilibrium_flux(t_idx, soln, gsp_inputs)

  subplot(132)
  plot_Bp(t_idx, soln, gsp_inputs)  

  subplot(133)
  plot_startup_metric(t_idx, soln, gsp_inputs)

  str = sprintf('Eq %d: Time=%.3f', t_idx, soln.t(t_idx));
  sgtitle(str, 'fontsize', 14, 'fontweight', 'bold')
  drawnow

end

function plot_startup_metric(t_idx, soln, gsp_inputs)

eq = soln.eqs{t_idx};
psizr = soln.signals.psizr;
tok = gsp_inputs.tok;
t = psizr.Time;

if t_idx == 1
  idx = [1 2];
else
  idx = [t_idx-1 t_idx];
end

% VACUUM
psidot = (soln.signals.psiapp.Data(idx(2),:) - soln.signals.psiapp.Data(idx(1),:)) /  (t(idx(2)) - t(idx(1)));
psiapp = soln.signals.psiapp.Data(t_idx,:);
psiapp = reshape(psiapp, tok.nz, tok.nr);
[~,~,Bp] = psizr2b(psiapp, tok);

% PLASMA
psidot = reshape(psidot, tok.nz, tok.nr);
Et = -psidot ./ (2*pi*tok.rgg);
Bt = eq.Btx;

startup_metric = Et.*Bt./Bp;
levels = quantile(startup_metric(~tok.outside_vessel_mask(:)), 0.05:0.15:0.95);
levels = round(levels, 2, 'significant');

cla
plot(tok.rl, tok.zl, 'k')
hold on
axis equal
axis([min(tok.rg) max(tok.rg) min(tok.zg) max(tok.zg)])
[C,h] = contour(tok.rg, tok.zg, startup_metric, levels, 'linewidth', 1); 
clabel(C,h,'labelspacing', 1000, 'fontsize', 14);
title('Startup metric Et*Bt/Bp [V/m]', 'fontweight', 'normal', 'fontsize', 14)


end


function plot_Bp(t_idx, soln, gsp_inputs)

  tok = gsp_inputs.tok;
  eq = soln.eqs{t_idx};
  
  cla
  plot(tok.rl, tok.zl, 'k')
  hold on
  axis equal
  axis([min(tok.rg) max(tok.rg) min(tok.zg) max(tok.zg)])

  Bp = sqrt(eq.Brx.^2 + eq.Bzx.^2) * 1e4;  % Gs
  levels = quantile(Bp(~tok.outside_vessel_mask(:)), 0.05:0.15:0.95);
  levels = round(levels, 2, 'significant');
  [C,h] = contour(tok.rg, tok.zg, Bp, levels, 'linewidth', 1); 
  clabel(C,h,'labelspacing', 1000, 'fontsize', 14);
  title('Bp [Gs]', 'fontweight', 'normal', 'fontsize', 14)
end



function plot_equilibrium_flux(t_idx, soln, gsp_inputs)

  tok = gsp_inputs.tok;
  shapes = gsp_inputs.shapes;
  eqs = soln.eqs;
  L = gsp_inputs.L;

  t = soln.t(t_idx);
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
  [~,~,~,~,~,~,rX,zX,FX] = meqpdom(eqs{t_idx}.Fx,eqs{t_idx}.Ip,1,L);
  plot(rX, zX, 'xk', 'linewidth', 2, 'markersize', 14)

  psi_max = max(eqs{t_idx}.Fx(~tok.outside_vessel_mask(:)));
  psi_min = min(eqs{t_idx}.Fx(~tok.outside_vessel_mask(:)));

  levels = linspace(psi_min, psi_max, 8);
  levels = [levels(:); FX(:)];
  contour(L.G.rx, L.G.zx, eqs{t_idx}.Fx, levels, 'color', [1 1 1] * 0.7, 'linewidth', 0.5)

  % plot equilibrium
  plot(tok.rl, tok.zl, 'k')
  if ~isempty(eqs{t_idx}.FB)
    contour(tok.rg, tok.zg, eqs{t_idx}.Fx, [1 1]*eqs{t_idx}.FB, 'r', 'linewidth', 1.5);  
  end

  % plot control points
  scatter(ref.cp_r, ref.cp_z, 20, 'b', 'filled')  
  title('Equilibrium', 'fontweight', 'normal', 'fontsize', 14)
  legend('Target', 'Actual', 'fontsize', 14)
  grid on
  box on

end
