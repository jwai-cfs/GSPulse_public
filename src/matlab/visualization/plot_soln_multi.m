function plot_soln_multi(soln1, soln2, gsp_inputs)
% =========================================================================
% Description: 
%  makes a series of overview plots of the GSPulse solution
%
% Inputs: 
%  soln1 - GSPulse solution struct
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

tlim = soln1.t([1 end]);
ccnames = tok.ccnames;
for i = 1:length(ccnames)
  coil_current_limits.(ccnames{i}).Time = tlim;
  coil_current_limits.(ccnames{i}).Data = [1 1]' * [settings.ic_min(i) settings.ic_max(i)];
  coil_voltage_limits.(ccnames{i}).Time = tlim;
  coil_voltage_limits.(ccnames{i}).Data = [1 1]' * [settings.vmin(i) settings.vmax(i)];
end

% Plot coil currents
axes(uitab(tg));
coil_currents1 = split_ts_signals(soln1.signals.ic, ccnames);
coil_currents2 = split_ts_signals(soln2.signals.ic, ccnames);
plot_structts(coil_currents1, ccnames, fig, 3, [], 'color', co{1});
plot_structts(coil_currents2, ccnames, fig, 3, [], 'color', co{2});
plot_structts(coil_current_limits, ccnames, fig, 3, [], '--r');
sgtitle('Coil Currents')

% Plot coil voltages
axes(uitab(tg));
coil_voltages1 = split_ts_signals(soln1.signals.voltage, ccnames);
coil_voltages2 = split_ts_signals(soln2.signals.voltage, ccnames);
plot_structts(coil_voltages1, ccnames, fig, 3, [], 'color', co{1});
plot_structts(coil_voltages2, ccnames, fig, 3, [], 'color', co{2});
plot_structts(coil_voltage_limits, ccnames, fig, 3, [], '--r');
sgtitle('Coil Voltages')

% Plot optimization signals
for i = 1:length(optimization_signals)
  sig = optimization_signals{i};
  if ismember(sig.name, cv.y_names)
    axes(uitab(tg))
    hold on
    plot(soln1.signals.(sig.name).Time, soln1.signals.(sig.name).Data, 'color', co{1})
    plot(soln2.signals.(sig.name).Time, soln2.signals.(sig.name).Data, 'color', co{2})
    plot(sig.target.Time, sig.target.Data, '--', 'color', co{3})
    title(strrep(sig.name, '_', ' '))
  end
end


%%
% Slider plot
slider_plot_soln_multi(soln1, soln2, gsp_inputs);
