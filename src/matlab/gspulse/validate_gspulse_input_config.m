function config = validate_gspulse_input_config(config)

% TODO: support defining tok directly with picard_algo='gspulse' instead of
% deriving from L
config.tok = L2tok(config.L);
assert(numel(config.initial_condition.x1) == config.tok.nc + config.tok.nv, ...
  'Number of coils (%d) and vessel (%d) modes not consistent with dimension of initial condition (%d).', ...
  config.tok.nc, config.tok.nv, numel(config.initial_condition.x1));

% validate combinations of the global settings
config.settings = validate_modify_settings(config.settings);

% validate shaping/dimensions of timeseries data
config.plasma_params = check_structts_dims(config.plasma_params);
config.shapes = check_structts_dims(config.shapes);
for i = 1:length(config.optimization_signals)
  config.optimization_signals{i} = check_structts_dims(config.optimization_signals{i});
end

% grab some timing information
time_data = parse_time_data(config.settings.interval_t);
t1 = time_data.interval(1).t1;
tN = time_data.interval(end).tN;
t0 = time_data.interval(1).t0;
tNm1 = time_data.interval(end).tNm1;

% validate timing of plasma params
config.plasma_params = validate_structts_timing(config.plasma_params, [t1 tN]);
config.plasma_params = remove_empty_fields(config.plasma_params);

% TODO: validate that optimization_signals contains exactly one signal with
% calc_type='voltage'

% TODO: validate that exactly one signal corresponds to 'coil_currents' and
% is named 'ic', and one signal for 'psibry'

% List which signals are being used (nonzero weight)
if config.settings.verbose > 1
  fprintf('\nSignals that have nonzero weights:\n')
  for i = 1:length(config.optimization_signals)
    sig = config.optimization_signals{i};
    if ~all(sig.wt.Data==0)
      fprintf('    %s.wt\n', sig.name);
    end
    if ~all(sig.dwt.Data==0)
      fprintf('    %s.dwt\n', sig.name);
    end
    if ~all(sig.d2wt.Data==0)
      fprintf('    %s.d2wt\n', sig.name);
    end
  end
  fprintf('\n')
end

% validate timing of shape targets
config.shapes = validate_structts_timing(config.shapes, [t1 tN]);

% TODO: out-of-bounds check that all shape targets lie within grid

% check consistency of initial condition, equality constraints
[~,~,~,targs,islocked] = assemble_weights_and_targets(config.optimization_signals);

for i = 1:time_data.n_intervals
  % check voltage constraints
  ilock = structts2vec(islocked, {'voltage'}, time_data.interval(i).t0Nm1);
  assert(all(ismember(ilock, [0 1])), 'Voltage constraint can only be within {0,1}. Interpolation produced non-logical value.')
  
  % check current constraints
  ilock = structts2vec(islocked, {'ic'}, time_data.interval(i).t1N);
  assert(all(ismember(ilock, [0 1])), 'Coil current constraint can only be within {0,1}. Interpolation produced non-logical value.')
end

% current constraints must be consistent with initial condition
if config.settings.verbose > 1
  fprintf('Checking constraint consistency coil currents\n'); 
end
k = 1:config.L.G.na;
vals = check_constraint_consistency(...
    config.initial_condition.x1(k), ...
    ~isnan(config.initial_condition.x1(k)), ...
    structts2vec(targs, {'ic'}, time_data.interval(1).t1), ...
    logical(structts2vec(islocked, {'ic'}, time_data.interval(1).t1)));
config.initial_condition.x1(k) = vals;

% voltage constraints must be consistent with initial condition
if config.settings.verbose > 1
  fprintf('Checking constraint consistency voltages\n\n'); 
end
config.initial_condition.u0 = check_constraint_consistency(...
    config.initial_condition.u0, ...
    ~isnan(config.initial_condition.u0), ...
    structts2vec(targs, {'voltage'}, time_data.interval(1).t0), ...
    logical(structts2vec(islocked, {'voltage'}, time_data.interval(1).t0)));
 
config.initial_condition.u1 = check_constraint_consistency(...
    config.initial_condition.u1, ...
    ~isnan(config.initial_condition.u1), ...
    structts2vec(targs, {'voltage'}, time_data.interval(1).t1), ...
    logical(structts2vec(islocked, {'voltage'}, time_data.interval(1).t1)));

% TODO: check that equality constraints do not lie outside current/voltage
%       limit (inequality constraints)

% validate optimization signals
for i = 1:length(config.optimization_signals)

  sig = config.optimization_signals{i};
  if strcmp(sig.calc_type, 'voltage')
    tspan = [t0 tNm1];
  else
    tspan = [t1 tN];
  end

  % validate timing
  if config.settings.verbose > 1
    fprintf('Checking timing signal "%s"\n', sig.name);
    if i == length(config.optimization_signals)
      fprintf('\n')
    end
  end
  config.optimization_signals{i} = validate_structts_timing(sig, tspan);

  % check shape dimensions
  switch sig.calc_type
    case 'flux_relative'
      assert(isequal( size(sig.r1.Data,2), size(sig.z1.Data,2), size(sig.targ.Data,2), ...
        size(sig.wt.Data,2), size(sig.dwt.Data,2), size(sig.d2wt.Data,2)), ...
        'Weight, target, and (r,z) dimensions are inconsistent for signal "%s"', sig.name)

      assert(isequal(1, size(sig.r2.Data,2), size(sig.z2.Data,2)))
    
    case 'flux_relative_multipoint'
      assert(isequal( size(sig.r.Data,2), size(sig.z.Data,2), size(sig.targ.Data,2), ...
        size(sig.wt.Data,2), size(sig.dwt.Data,2), size(sig.d2wt.Data,2)), ...
        'Weight, target, and (r,z) dimensions are inconsistent for signal "%s"', sig.name)
      
      assert(isequal( size(sig.r_multiref.Data,2), size(sig.z_multiref.Data,2), size(sig.wt_multiref.Data,2)), ...
        'Multiref weight, and (r,z) dimensions are inconsistent for signal "%s"', sig.name)
      
      assert(all(sum(sig.wt_multiref.Data, 2)==1), ...
        'Multiref weights must sum to 1 at all times')


    case {'flux_absolute', 'field_absolute_vertical', 'field_absolute_radial', ...
        'vacuum_field_absolute_vertical', 'vacuum_field_absolute_radial'}
      
      assert(isequal( size(sig.r.Data,2), size(sig.z.Data,2), size(sig.targ.Data,2), ...
        size(sig.wt.Data,2), size(sig.dwt.Data,2), size(sig.d2wt.Data,2)), ...
        'Weight, target, and (r,z) dimensions are inconsistent for signal "%s"', sig.name)

    case {'flux_absolute_avg', 'vacuum_flux_absolute_avg'}
      assert(isequal( 1, size(sig.targ.Data,2), ...
        size(sig.wt.Data,2), size(sig.dwt.Data,2), size(sig.d2wt.Data,2)), ...
        'Weight, target, and (r,z) dimensions are inconsistent for signal "%s"', sig.name)
    
    case 'coil_current_combinations'
      assert(isequal(size(sig.coil_combinations_matrix,1), size(sig.targ.Data,2), ...
        size(sig.wt.Data,2), size(sig.dwt.Data,2), size(sig.d2wt.Data,2)), ...
        'Weight, target, and combo_matrix dimensions are inconsistent for signal "%s"', sig.name)

      assert(isequal(config.tok.nc, size(sig.coil_combinations_matrix,2)), ['Number of columns of ' ...
        'combo matrix must match number of coils, for signal "%s"'], sig.name)

    case 'coil_currents'
      assert(isequal( config.tok.nc, size(sig.targ.Data,2), ...
        size(sig.wt.Data,2), size(sig.dwt.Data,2), size(sig.d2wt.Data,2)), ...
        'Weight, and target dimensions must be consistent with # coils, for signal "%s"', sig.name)
    
    case 'voltage'
      assert(isequal( config.tok.nc, ...
        size(sig.wt.Data,2), size(sig.dwt.Data,2), size(sig.d2wt.Data,2)), ...
        'Dimensions of weights are inconsistent with # coils for signal "%s"', sig.name)
    
    otherwise
      error('Unknown calculation type %s', sig.calc_type)
  end
end

end

%% local functions
function [vals, ilock] = check_constraint_consistency(vals1,ilock1,vals2,ilock2)
  bothlock = ilock1 & ilock2;  
  assert( isequal(vals1(bothlock), vals2(bothlock)), ...
    'Inconsistency between the initial condition constraint and per-coil waveform constraints.')
  
  vals = nan(size(vals1)); 
  vals(ilock1) = vals1(ilock1); 
  vals(ilock2) = vals1(ilock2);   
  ilock = ilock1 | ilock2;  
end