function fig = slider_plot_soln_multi(soln1, soln2, gsp_inputs)

% create figure
fig = figure;
fig.Position = [1044 491 444 705];
ax = axes(fig);
ax.Box = 'on';
ax.Position = [0.13 0.2 0.775 0.7];

% slider button
s = uicontrol(fig, 'style', 'slider');
s.Units = 'normalized';
s.Position = [0.15 0.05 0.6 0.05];
s.Min = 1;
s.Max = length(soln1.t);
s.Value = 1;
s.SliderStep = [1 1] / (s.Max - s.Min);
s.Callback = {@sliderCallback, soln1, soln2, gsp_inputs};

update_plots(soln1.t(1), soln1, soln2, gsp_inputs)

% slider callback
function sliderCallback(src, ~, soln1, soln2, gsp_inputs)
  t_idx = src.Value;
  t_idx = round(t_idx);
  t = soln1.t(t_idx);
  update_plots(t, soln1, soln2, gsp_inputs)
end

end

% update plots
function update_plots(t, soln1, soln2, gsp_inputs)
    
  cla
  co = mat2cell(colororder, ones(7,1), 3);
  tok = gsp_inputs.tok;
  plot(tok.rl, tok.zl, 'k')
  axis equal
  axis([min(tok.rg) max(tok.rg) min(tok.zg) max(tok.zg)])
  hold on

  plot_equilibrium_flux(t, soln1, gsp_inputs, co{1})
  plot_equilibrium_flux(t, soln2, gsp_inputs, co{2})
  str = sprintf('Time=%.3f', t);
  sgtitle(str, 'fontsize', 14, 'fontweight', 'bold')
  drawnow

end


function plot_equilibrium_flux(t, soln, gsp_inputs, color)

  [~,t_idx] = min(abs(t-soln.t));
  tok = gsp_inputs.tok;
  shapes = gsp_inputs.shapes;
  eqs = soln.eqs;
  L = gsp_inputs.L;
 
  shape_fields = fieldnames(shapes);
  ref = structts2struct(shapes, shape_fields, t, 0);
  
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
  [~,~,~,~,~,~,rX,zX,FX] = meqpdom(eqs{t_idx}.psizr,eqs{t_idx}.Ip,1,L);
  plot(rX, zX, 'xk', 'linewidth', 2, 'markersize', 14)

  psi_max = max(eqs{t_idx}.psizr(~tok.outside_vessel_mask(:)));
  psi_min = min(eqs{t_idx}.psizr(~tok.outside_vessel_mask(:)));
  
  % plot equilibrium  
  if ~isempty(eqs{t_idx}.FB)
    contour(tok.rg, tok.zg, eqs{t_idx}.psizr, [1 1]*eqs{t_idx}.FB, 'color', color, 'linewidth', 1.5);  
  end

  % plot control points
  scatter(ref.cp_r, ref.cp_z, 20, 'b', 'filled')  
  title('Equilibrium', 'fontweight', 'normal', 'fontsize', 14)  
  grid on
  box on

end
