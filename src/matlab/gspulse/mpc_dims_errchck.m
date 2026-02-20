function mpc_dims_errchck(settings, shapes, weights, targs, Cmats, cv)
% =========================================================================
% Description: 
%  check for dimensional and timing consistency betweeen the shapes,
%  weights, and targs
%
% Inputs: 
%
%  settings - settings struct, see help _define_settings.m
%  shapes   - shapes struct, see help _define_shapes.m
%  weights  - weights struct, see help _define_weights.m
%  targs    - targets struct, see help _define_targets.m
%  Cmats    - output model linearizations, see output_model.m
%  cv       - data indices, see define_data_indices.m 
% 
% Outputs: 
%  None - prints warnings about dimensional and timing consistency
%
% =========================================================================
for fd_ = settings.fds2control(:)'
  fd = fd_{:};
 
  % check dimensions for targs and weights
  ntargs = size(targs.(fd).Data, 2); 
  nwts = size(weights.wts.(fd).Data, 2);
  ndwts = size(weights.dwts.(fd).Data, 2);
  nd2wts = size(weights.d2wts.(fd).Data, 2);  
  if ntargs ~= nwts    
    warning(['Dimension inconsistency targs.%s vs weights.wts.%s, ' ...
      'expected %d got %d'], fd, fd, ntargs, nwts);
  end  
  if ntargs ~= ndwts    
    warning(['Dimension inconsistency targs.%s vs weights.dwts.%s, ' ...
      'expected %d got %d'], fd, fd, ntargs, ndwts);
  end  
  if ntargs ~= nd2wts    
    warning(['Dimension inconsistency targs.%s vs weights.d2wts.%s, ' ...
      'expected %d got %d'], fd, fd, ntargs, nd2wts);
  end
  
  % check response model for nans
  i = cv.iy.(fd);
  found_nan = 0;
  for j = 1:length(Cmats)    
    if any(isnan(Cmats{j}(i,:)))
      found_nan = 1;
    end
  end
  if found_nan
    warning('Output response model (Cmats) contains nan for "%s" response', fd);
  end

  % % check that data exists in time window
  % t0 = settings.interval_t{1}(1);
  % tf = settings.interval_t{end}(end);
  % if weights.wts.(fd).Time(1) > t0 || weights.wts.(fd).Time(end) < tf
  %   warning(['weights.wts.%s not specified for entire time interval ' ...
  %     '%.3f to %.3f, endpoints will be extrapolated'], fd, t0, tf);
  % end
  % if weights.dwts.(fd).Time(1) > t0 || weights.dwts.(fd).Time(end) < tf
  %   warning(['weights.dwts.%s not specified for entire time interval ' ...
  %     '%.3f to %.3f, endpoints will be extrapolated'], fd, t0, tf);
  % end
  % if weights.d2wts.(fd).Time(1) > t0 || weights.d2wts.(fd).Time(end) < tf
  %   warning(['weights.d2wts.%s not specified for entire time interval ' ...
  %     '%.3f to %.3f, endpoints will be extrapolated'], fd, t0, tf);
  % end
  % if targs.(fd).Time(1) > t0 || targs.(fd).Time(end) < tf
  %   warning(['targs.%s not specified for entire time interval ' ...
  %     '%.3f to %.3f, endpoints will be extrapolated'], fd, t0, tf);
  % end

  % check if targs contains nans
  found_nan = 0;
  for j = 1:length(settings.interval_t)
    t = settings.interval_t{j};
    data = structts2vec(targs, {fd}, t);
    if any(isnan(data))
      found_nan = 1;
    end
  end
  if found_nan
    warning('targs.%s contains nan.', fd)
  end

  % check consistency of shapes and targs
  chk1 = {'rb', 'zb'};
  chk2 = {'diff_psicp_psitouch', 'diff_psicp_psix1', 'diff_psicp_psix2'};
  chk2 = intersect(chk2, settings.fds2control(:)');
  for a = chk1
    for b = chk2
      na = size(shapes.(a{:}).Data, 2);
      nb = size(targs.(b{:}).Data, 2);
      if  na ~= nb
        warning('Dimension of shapes.%s (%d) does not match targs.%s (%d)', a{:}, na, b{:}, nb);
      end
    end
  end

  chk1 = {'rstrike', 'zstrike'};
  chk2 = {'diff_psisp_psitouch', 'diff_psisp_psix1', 'diff_psisp_psix2'};
  chk2 = intersect(chk2, settings.fds2control);
  for a = chk1
    for b = chk2
      na = size(shapes.(a{:}).Data, 2);
      nb = size(targs.(b{:}).Data, 2);
      if  na ~= nb
        warning('Dimension of shapes.%s (%d) does not match targs.%s (%d)', a{:}, na, b{:}, nb);
      end
    end
  end
end