function settings = validate_modify_settings(settings)

% check string inputs
% TODO: make these inputs enums and remove check
settings = enum_check(settings, 'qpsolver', {'quadprog', 'scs', 'python_cvxopt'});
settings = enum_check(settings, 'picard_algo', {'fbt'});  % TODO: support, picard_algo = 'gspulse' 
settings = enum_check(settings, 'specify_psibry_mode', {'ejima', 'direct'});
 
tmp_ = dir([getenv('GSROOT') '/tokamaks']);
toks = {tmp_(:).name};
toks = toks(~startsWith(toks, '.'));
settings = enum_check(settings, 'tokamak', toks);

% octave compatibility
is_octave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if is_octave  
  oct = struct;
  oct.do_final_boundary_trace.allow = {0}; 
  oct.do_final_boundary_trace.default = 0; 
  oct.plotlevel.allow = {0};
  oct.plotlevel.default = 0; 
  oct.qpsolver.allow = {'python_cvxopt', 'scs'};
  oct.qpsolver.default = 'python_cvxopt';
  
  fds = fields(oct);
  for i = 1:length(fds)
    fd = fds{i};
    if ~iscellmember(settings.(fd), oct.(fd).allow)         
      warning('settings.%s="%s" not supported in octave (allowed={%s}), setting to octave default "%s".\n', ...
        fd, tostring(settings.(fd)), tostring(oct.(fd).allow), tostring(oct.(fd).default))
      settings.(fd) = oct.(fd).default;
    end
  end
end

% silent overrides
overrides.dr_OMP = 0.001; 
overrides.direction = 1;
overrides.npts = 101; 
fds = fields(overrides);
for i = 1:length(fds)
    fd = fds{i};
  if ~isfield(settings, fd)
    settings.(fd) = overrides.(fd);
  end
end
 
% very specific settings combinations
% --------------------------------------

% additional rec on psibry
if strcmpi(settings.specify_psibry_mode, 'ejima') && length(settings.interval_t) > 1 
  warning(['Recommended to specify target psibry directly, computing target ' ...
    'psibry automatically when using multiple time intervals is not well ' ...
    'supported.']);
end

% consistency of timing and splining
if settings.use_spline_basis 
  for i = 1:length(settings.interval_t)
    if length(settings.interval_t{i}) < 3
      warning(['Using a spline basis for the solve is not supported when ' ...
        'there are fewer than 3 time points in a time interval. Modifying the ' ...
        'value such that settings.use_spline_basis = 0'])
      settings.use_spline_basis = 0;
    end
  end
end

% consistency of timing and using multiple intervals for solve
if length(settings.interval_t) > 1
  for i = 1:length(settings.interval_t)
    if length(settings.interval_t{i}) < 3
      warning(['Each time interval must contain 3 or more time points, when ' ...
        'using multiple time intervals.'])      
        settings.use_spline_basis = 0;
    end
  end
end

% strike point calc only supports sparc currently
if settings.calc_strike_pts && ~strcmpi(settings.tokamak, 'sparc')
  warning("calc_strike_pts=true only supported for tokamak='sparc', setting to false.")
  settings.calc_strike_pts = 0;
end

assert(settings.Iy_relax_factor <= 1 && settings.Iy_relax_factor >= 0, ...
  "Iy_relax_factor must be within {0,1}, received %f\n", settings.Iy_relax_factor)

% unused fields
used_fds = {'pulse_id', 'tokamak', 'picard_algo', 'specify_psibry_mode', ... 
  'qpsolver', 'qpsolver_tol', 'spline_basis_ratio', 'plotlevel', 'niter', ...
  'verbose', 'do_final_boundary_trace', 'use_spline_basis', 'calc_strike_pts', ...
  'inject_model_offset', 'calc_post_proc_ext', 'ncoils', 'ncombos', ...
  'interval_t', 'ic_max', 'ic_min', 'vmax', 'vmin', 'pulse_id_str', ...
  'calc_post_prc_ext', 'tol_ic_diff', 'boundary_post_run', 'dr_OMP', ...
  'direction', 'npts', "fbt_Ip_threshold", "Iy_smooth_dt", "Iy_relax_factor"};
fds = fields(settings);
for i = 1:length(fds)
  fd = fds{i};
  if ~iscellmember(fd, used_fds)
    warning('Settings contains the unused field "%s"', fd);
  end
end

end

%% local functions
% --------------------------------
function settings = enum_check(settings, nam, options)
  assert( any(strcmpi(settings.(nam), options)), ...
    "Unsupported value settings.%s = %s, must be one of {%s}.\n", ...
    nam, tostring(settings.(nam)), tostring(options))
  settings.(nam) = lower(settings.(nam));
end

function ismember = iscellmember(m,c)
  ismember = false;
  for i = 1:numel(c)
    if isequal(c{i}, m)
      ismember=true; 
      break
    end
  end
end

function str = tostring(a)
  if ischar(a)
    str = a;
  elseif isnumeric(a)
    str = num2str(a);
  elseif iscell(a)  
    str = '';
    for i = 1:numel(a)
      str = [tostring(a{i}), ',' str];
    end
    str(end) = [];
  end
end

