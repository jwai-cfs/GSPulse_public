function plot_soln(soln, gsp_inputs)
% =========================================================================
% Description: 
%  makes a series of overview plots of the GSPulse solution
%
% Inputs: 
%  soln - GSPulse solution struct
%  gsp_inputs - the GSPulse inputs read from config files, also see
%               run_pulse.m
%
% Outputs: 
%  None, makes a bunch of plots
%
% =========================================================================

settings = gsp_inputs.settings;
tok = gsp_inputs.tok;
shapes = gsp_inputs.shapes;
optimization_signals = gsp_inputs.optimization_signals;
co = mat2cell(colororder, ones(7,1), 3);
cv = define_data_indices(optimization_signals,tok);

fig = figure;
fig.Position = [342 427 560 420];
tg = uitabgroup;
set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')

tlim = soln.t([1 end]);
ccnames = tok.ccnames;
for i = 1:length(ccnames)
  coil_current_limits.(ccnames{i}).Time = tlim;
  coil_current_limits.(ccnames{i}).Data = [1 1]' * [settings.ic_min(i) settings.ic_max(i)];
  coil_voltage_limits.(ccnames{i}).Time = tlim;
  coil_voltage_limits.(ccnames{i}).Data = [1 1]' * [settings.vmin(i) settings.vmax(i)];
end


if numel(soln.stage) > 1

  % Plot overlapping stage coil currents
  axes(uitab(tg));
  for i = 1:numel(soln.stage)
    coil_currents = split_ts_signals(soln.stage{i}.mpcsoln.ic, ccnames);
    if mod(i,2)
      plot_structts(coil_currents, ccnames, fig, 3, [], '-ob', 'markersize', 5);
    else
      plot_structts(coil_currents, ccnames, fig, 3, [], '-or', 'markersize', 3, 'markerfacecolor', 'r');
    end
    plot_structts(coil_current_limits, ccnames, fig, 3, [], '--r');
    set(gca, 'xlim', soln.t([1 end]))
    sgtitle('Coil currents for each optimization interval')
  end
  
  % Plot overlapping stage voltages
  axes(uitab(tg));
  for i = 1:numel(soln.stage)
    coil_voltages = split_ts_signals(soln.stage{i}.mpcsoln.voltage, ccnames);
    if mod(i,2)
      plot_structts(coil_voltages, ccnames, fig, 3, [], '-ob', 'markersize', 5);
    else
      plot_structts(coil_voltages, ccnames, fig, 3, [], '-or', 'markersize', 3, 'markerfacecolor', 'r');
    end
    plot_structts(coil_voltage_limits, ccnames, fig, 3, [], '--r');
    set(gca, 'xlim', soln.t([1 end]))
    sgtitle('Coil voltages for each optimization interval')
  end
end

% Plot coil currents
axes(uitab(tg));
coil_currents = split_ts_signals(soln.signals.ic, ccnames);
plot_structts(coil_currents, ccnames, fig, 3);
plot_structts(coil_current_limits, ccnames, fig, 3, [], '--r');
sgtitle('Coil Currents')

% Plot coil voltages
axes(uitab(tg));
coil_voltages = split_ts_signals(soln.signals.voltage, ccnames);
plot_structts(coil_voltages, ccnames, fig, 3);
plot_structts(coil_voltage_limits, ccnames, fig, 3, [], '--r');
sgtitle('Coil Voltages')

% Plot optimization signals
for i = 1:length(optimization_signals)
  sig = optimization_signals{i};
  if ismember(sig.name, cv.y_names)
    axes(uitab(tg))
    hold on
    plot(soln.signals.(sig.name).Time, soln.signals.(sig.name).Data, 'color', co{1})
    plot(sig.targ.Time, sig.targ.Data, '--', 'color', co{2})
    title(strrep(sig.name, '_', ' '))
  end
end


% Plot startup signal
axes(uitab(tg))

% on-midplane Et, Bt, Bp, Br, Bz
r = [1.6 1.7 1.8 1.9]'; 
z = zeros(size(r));
[psi, psi_r, psi_z] = multigrid2pt(tok.rg, tok.zg, soln.signals.psiapp.Data', r, z);
Bz = diag(1./(2*pi*r)) * psi_r;
Br = -diag(1./(2*pi*r)) * psi_z;
Bp = sqrt(Br.^2 + Bz.^2);
Bt = diag(soln.eqs{1}.rBt ./ r);
psidot = diff(psi,1,2) * diag(1./diff(soln.t));
psidot = psidot(:, [1 1:end]);
Et = -diag(1 ./ (2*pi*r)) * psidot;
su_metric = Bt*Et ./ Bp;
labels_on_midplane = {'r=1.6m,z=0', 'r=1.7m,z=0', 'r=1.8m,z=0', 'r=1.9m,z=0'};

% off-midplane Br
r = [1.6 1.7 1.8]'; 
z = ones(size(r)) * 0.6;
[~,~,psi_z] = multigrid2pt(tok.rg, tok.zg, soln.signals.psiapp.Data', r, z);
Br = -diag(1./(2*pi*r)) * psi_z;
labels_off_midplane = {'r=1.6m,z=0.6m', 'r=1.7m,z=0.6m', 'r=1.8m,z=0.6m'};

% plots
subplot(221)
plot(soln.t, Et)
title('Vacuum Et [V/m]', 'fontweight', 'normal')
legend(labels_on_midplane)
subplot(222)
plot(soln.t, Bz)
title('Vacuum Bz [T]', 'fontweight', 'normal')
legend(labels_on_midplane)
subplot(223)
plot(soln.t, su_metric)
ylim([-100 2000])
title('Vacuum Et*Bt/Bp [V/m]', 'fontweight', 'normal')
legend(labels_on_midplane)
subplot(224)
plot(soln.t, Br)
title('Vacuum Br [T]', 'fontweight', 'normal')
legend(labels_off_midplane)

%%
% Slider plot
slider_plot_plasma_startup(soln, gsp_inputs);
