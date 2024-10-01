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

if settings.plotlevel > 0

  fig = figure;
  tg = uitabgroup;
  set(groot,'defaultAxesXGrid','on')
  set(groot,'defaultAxesYGrid','on')

  tlim = soln.t([1 end]);
  ctags = tok.ccnames;
  vtags = strcat(ctags, '_V');
  for i = 1:length(ctags)
    limits.(ctags{i}).Time = tlim;
    limits.(ctags{i}).Data = [1 1]' * [settings.ic_min(i) settings.ic_max(i)];
    limits.(vtags{i}).Time = tlim;
    limits.(vtags{i}).Data = [1 1]' * [settings.vmin(i) settings.vmax(i)];
  end

  % plot overlapping stage currents
  if 0
    axes(uitab(tg));
    plot_structts(soln.stage{1}.mpcsoln, tok.ccnames(1:2:end), fig, 3, [], '-r.', 'markersize', 10);  % plots individual coil currents
    plot_structts(soln.stage{2}.mpcsoln, tok.ccnames(1:2:end), fig, 3, [], '-bo');
    sgtitle('Coil currents')
    set(gcf, 'Position', [819 436 554 428])
    drawnow
  end
  
  % plot overlapping stage voltages
  if 0
    axes(uitab(tg));
    targs = soln.targs;
    cv = define_data_indices(settings, targs, tok);
    coils = fieldnames(cv.iu);
    coils = coils(1:2:end);
    plot_structts(soln.stage{1}.mpcsoln, coils, fig, 3, [], '-r.', 'markersize', 10);  % plots individual coil currents
    plot_structts(soln.stage{2}.mpcsoln, coils, fig, 3, [], '-bo');
    sgtitle('Coil voltages')
    set(gcf, 'Position', [254 441 564 425])
    drawnow
  end
  
  % plot coil currents
  if 1
    axes(uitab(tg));
    plot_structts(soln.signals, ctags, fig, 3);
    plot_structts(limits                 , ctags, fig, 3, [], '--r');
    sgtitle('Coil Currents')
  end
    
  % plot coil voltages
  if 1
    axes(uitab(tg));
    plot_structts(soln.signals, vtags, fig, 3);
    plot_structts(limits                 , vtags, fig, 3, [], '--r');
    sgtitle('Coil voltages')
  end
  
  % psibry
  if 1
    axes(uitab(tg));
    x = merge_structs(soln.targs, soln.plasma_params);
    plot_structts(soln.signals, {'psibry'}, fig, 2, 2, 'b', 'linewidth', 2);
    hold on
    fds2plot = {'psibry', 'Wk', 'Ip'};
    plot_structts(x, fds2plot, fig, 2, 2, '--r');    
  end

  % shape slider plot
  if 1
    shape_slider_plot(soln.t, shapes, soln.eqs, tok);  
  end

end


 