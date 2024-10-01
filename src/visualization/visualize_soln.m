if settings.plotlevel > 0

  fig = figure;
  tg = uitabgroup;
  set(groot,'defaultAxesXGrid','on')
  set(groot,'defaultAxesYGrid','on')

  tlim = soln.globalsoln.t([1 end]);
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
    try
    axes(uitab(tg));
    plot_structts(soln.stagesoln{1}.mpcsoln, tok.ccnames(1:2:end), fig, 3, [], '-r.', 'markersize', 10);  % plots individual coil currents
    plot_structts(soln.stagesoln{2}.mpcsoln, tok.ccnames(1:2:end), fig, 3, [], '-bo');
    sgtitle('Coil currents')
    xlim([0 10])
    set(gcf, 'Position', [819 436 554 428])
    drawnow
    catch
    end
  end
  
  % plot overlapping stage voltages
  if 0
    try
    axes(uitab(tg));
    cv = define_data_indices(settings, targs, tok);
    coils = fieldnames(cv.iu);
    coils = coils(1:2:end);
    plot_structts(soln.stagesoln{1}.mpcsoln, coils, fig, 3, [], '-r.', 'markersize', 10);  % plots individual coil currents
    plot_structts(soln.stagesoln{2}.mpcsoln, coils, fig, 3, [], '-bo');
    sgtitle('Coil voltages')
    xlim([0 10])
    set(gcf, 'Position', [254 441 564 425])
    drawnow
    catch
    end
  end
  
  % plot coil currents
  if 1
    axes(uitab(tg));
    plot_structts(soln.globalsoln.mpcsoln, ctags, fig, 3);
    plot_structts(limits                 , ctags, fig, 3, [], '--r');
    sgtitle('Coil Currents')
  end
    
  % plot coil voltages
  if 1
    axes(uitab(tg));
    plot_structts(soln.globalsoln.mpcsoln, vtags, fig, 3);
    plot_structts(limits                 , vtags, fig, 3, [], '--r');
    sgtitle('Coil voltages')
  end
  
  % psibry
  if 1
    axes(uitab(tg));
    x = merge_structs(soln.globalsoln.targs, soln.plasma_params);
    plot_structts(soln.globalsoln.mpcsoln, {'psibry'}, fig, 2, 2);
    plot_structts(x, {'psibry', 'Wk', 'Ip'}, fig, 2, 2, '--r');    
  end

  % shape slider plot
  if 1
    shape_slider_plot(soln.globalsoln.t, shapes, soln.globalsoln.eqs, tok);  
  end

end


 