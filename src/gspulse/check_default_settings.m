function modified_settings = check_default_settings(settings, tokamak)
% =========================================================================
% Description: 
% check that settings defined by the user have been specified correctly,
% and assign default values to the optional fields
%
% 
% Inputs: 
% settings - settings struct, see help _define_general_settings.m
%
% Outputs: 
% modified_settings - settings that have been checked/formatted/assigned
%                     default values
%
% =========================================================================

settings.tokamak = tokamak;
modified_settings = settings;

% required settings - must be provided by user
required = {'interval_t', 'fds2control', 'vmax', 'vmin', 'ic_max', ...
            'ic_min', 'tokamak'};
for x = required(:)'
  if ~isfield(settings, x{:})
    warning('settings does not contain required field "%s"', x{:});
  end
end

% default settings - will be provided if not supplied by user
d.plotlevel = 1;                 
d.verbose = 1;                   
d.dosave = false;                
d.niter = 10;                    
d.use_spline_basis = 0;
d.spline_basis_ratio = 0.7;
d.nvessmodes = 40;
d.enforce_voltage_limits = 1;
d.enforce_current_limits = 1;
d.picard_algo = 'fbt';
d.do_final_boundary_trace = 1;
d.inject_model_offset = 0;
d.specify_psibry_mode = 'direct';
d.qpsolver = 'quadprog'; % by default use MATLAB's quadprog
d.calc_strike_pts = 0;

% overwrite missing fields with defaults
defaults = fieldnames(d);
for x = defaults(:)'
  if ~isfield(settings, x{:})
    warning(['settings does not specify "%s", setting to default value ' ...
      'of "%s".'], x{:}, num2str(d.(x{:})));
    modified_settings.(x{:}) = d.(x{:});
  end  
end

settings = modified_settings;

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
      modified_settings.use_spline_basis = 0;
    end
  end
end

% consistency of timing and using multiple intervals for solve
if length(settings.interval_t) > 1
  for i = 1:length(settings.interval_t)
    if length(settings.interval_t{i}) < 3
      warning(['Each time interval must contain 3 or more time points, when ' ...
        'using multiple time intervals.'])
      modified_settings.use_spline_basis = 0;
    end
  end
end

% check that no additional fields are specified
settings_fds = [required'; fieldnames(d)];
check_fds = fieldnames(settings);
for i = 1:length(check_fds)
  check_fd = check_fds{i};
  if ~ismember(check_fd, settings_fds)
    warning('Settings contains the unused field "%s"', check_fd);
  end
end

% strike point calc only supports sparc currently
if settings.calc_strike_pts && ~strcmpi(settings.tokamak, 'sparc')
  warning("calc_strike_pts=true only supported for tokamak='sparc', setting to false.")
  settings.calc_strike_pts = 0;
end

