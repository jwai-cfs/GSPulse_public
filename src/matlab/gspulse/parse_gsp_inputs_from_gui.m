function gsp_inputs = parse_gsp_inputs_from_gui(gui_data_filepath, settings_overrides)
% ===========================================
% Build gsp_inputs from GUI-specified data
% ===========================================
if ~exist('settings_overrides', 'var'), settings_overrides = struct; end

% read GUI inputs
gui_data = jsondecode(fileread(gui_data_filepath));
settings = gui_data.settings; 
fds = fields(settings_overrides);
for i = 1:length(fds)
  settings.(fds{i}) = settings_overrides.(fds{i});
end
gsp_inputs = struct;

% settings
% -----------

% parse time interval
settings.interval_t = {};
for i = 1:4
  fd = ['interval_t' num2str(i)];
  t = settings.(fd);
  if ~isempty(t) || any(isnan(t))
    settings.interval_t{i} = t;
  end
  settings = rmfield(settings, fd);
end
settings.interval_t = settings.interval_t(:);
 
% read current/voltage limits
coils = fields(gui_data.ic_max);
settings.ic_max = struct2vec(gui_data.ic_max, coils) * 1e3;
settings.ic_min = struct2vec(gui_data.ic_min, coils) * 1e3;
settings.vmax = struct2vec(gui_data.vmax, coils) * 1e3;
settings.vmin = struct2vec(gui_data.vmin, coils) * 1e3;

% pulse id
id_str = settings.pulse_id;
id = str2double(id_str);
settings.pulse_id_str = id_str;
if isempty(id) || isnan(id)
  id = 0;
end
settings.pulse_id = id;

% misc cleanup
settings = recursive_bools_to_doubles(settings);
settings = rmfield(settings, 'initial_vess_currents');
gsp_inputs.settings = validate_modify_settings(settings);


% plasma_params
% -----------------
plasma_params = struct;

% grab some timing information
time_data = parse_time_data(gsp_inputs.settings.interval_t);
t1 = time_data.interval(1).t1;
tN = time_data.interval(end).tN;
t0 = time_data.interval(1).t0;
tNm1 = time_data.interval(end).tNm1;

% plasma scalar waveforms
for var_ = parse_fbtagcon_str(gui_data.fbt_params.fbtagcon)
  var = var_{:};   % e.g. var = 'Ip'
  plasma_params.(var).Time = gui_data.plasma_scalars.([var '_time']);
  plasma_params.(var).Data = gui_data.plasma_scalars.([var '_data']);
  if strcmp(var, 'ag')
    plasma_params.ag.Data = plasma_params.ag.Data * [1 1 1]; % ag formatting
  end
end
plasma_params.Ip.Data = plasma_params.Ip.Data*1e6;  % MA-->A
plasma_params.Wk.Data = plasma_params.Wk.Data*1e6;  % MJ-->J
plasma_params.rBt.Time = [t1; tN];
plasma_params.rBt.Data = [1;1] * gui_data.fbt_params.r0 * gui_data.fbt_params.b0;

% read profiles
npsin = 101;
plasma_params.pprime = interpolate_profile_to_psin(gui_data.pprime, npsin);
plasma_params.ttprime1 = interpolate_profile_to_psin(gui_data.ffprime1, npsin);
plasma_params.ttprime2 = interpolate_profile_to_psin(gui_data.ffprime2, npsin);

% integrate pprime for pressure
plasma_params.press.Time = plasma_params.pprime.Time;
for i = 1:length(plasma_params.press.Time)  
  plasma_params.press.Data(i,:) = bfprmex(plasma_params.pprime.Data(i,:)' * ones(1,3)) * [1 0 0]';    
end

% validate timing of plasma params
gsp_inputs.plasma_params = validate_structts_timing(plasma_params, [t1 tN]);


% shapes
% ------------
shapes = struct;

% get shape field names
fds = fieldnames(gui_data.shape_array_data); 
i = ~startsWith(fds, 'cp_') | (strcmp(fds, 'cp_r') | strcmp(fds, 'cp_z'));
old_shape_fdnames = fds(i);

% rename several variables
new_shape_fdnames = old_shape_fdnames;
new_shape_fdnames = strrep(new_shape_fdnames, 'R_strike_pt', 'rstrike');
new_shape_fdnames = strrep(new_shape_fdnames, 'Z_strike_pt', 'zstrike');
new_shape_fdnames = strrep(new_shape_fdnames, 'R_x_pt', 'rx');
new_shape_fdnames = strrep(new_shape_fdnames, 'Z_x_pt', 'zx');
new_shape_fdnames = strrep(new_shape_fdnames, 'R_control_pt_ref', 'r_control_pt_ref');
new_shape_fdnames = strrep(new_shape_fdnames, 'Z_control_pt_ref', 'z_control_pt_ref');

% read data
for i = 1:length(old_shape_fdnames)
  old_fd = old_shape_fdnames{i};
  new_fd = new_shape_fdnames{i};
  shapes.(new_fd).Time = gui_data.shape_array_time;
  shapes.(new_fd).Data = gui_data.shape_array_data.(old_fd)';
end
gsp_inputs.shapes = shapes;

% initial condition
% --------------------

% initial vessel currents
if isscalar(gui_data.settings.initial_vess_currents)
  iv1 = ones(gui_data.fbt_params.nu,1) * gui_data.settings.initial_vess_currents;
else
  iv1 = gui_data.settings.initial_vess_currents(:);
end

% initial coil currents
ic1 = recursive_replace_vals(gui_data.ic_init, [], nan);
ic1 = struct2vec(ic1, coils);

% initial voltages
u0 = recursive_replace_vals(gui_data.vinit, [], nan);
u0 = struct2vec(u0, coils);

% format
gsp_inputs.initial_condition = struct;
gsp_inputs.initial_condition.x1 = [ic1; iv1];
gsp_inputs.initial_condition.u0 = u0;
gsp_inputs.initial_condition.u1 = nan(size(u0));


% L geometry obj
% -----------------

% parse P overrides
P = gui_data.fbt_params;
switch P.bfct
  case 'bfabmex'
    P.bfp = [1 2];
  case 'bf3imex'
    % do nothing
  otherwise
    error('Unsupported value for "bfct"')
end
P.bfct = str2func(P.bfct);
P.fbtagcon = parse_fbtagcon_str(P.fbtagcon);
P.shot = id;

% load L
Lfun = str2func(sprintf('load_%s_L_fbt', settings.tokamak));
L = Lfun(P); 
L.G = meqg(L.G, L.P,'Maa','Mau','Mav','Muu','Ra','Ru');
gsp_inputs.L = L;
gsp_inputs.tok = L2tok(gsp_inputs.L);


% optimization signals 
% ---------------------- 
optimization_signals = {};

% voltage
if gui_data.voltage.use_in_optimization
  sig = struct; 
  sig.name = 'voltage';
  sig.calc_type = 'voltage';
  sig = parse_add_per_coil_waveforms(sig, gui_data.voltage, coils);
  % extrapolate to include t=t0 (before the first equilibrium solve)
  for fd_ = {'wt', 'dwt', 'd2wt', 'targ', 'constrain'}
    fd = fd_{:};
    t = sort([t0; sig.(fd).Time]);
    sig.(fd).Data = interp1hold(sig.(fd).Time, sig.(fd).Data, t);
    sig.(fd).Time = t;
  end
  sig.targ.Data = sig.targ.Data * 1e3; % kV --> V
  % by default, scale voltage weight by coil inductance
  scale = sparse(diag(diag(gsp_inputs.L.G.Maa))) .^2;
  sig.wt.Data = sig.wt.Data * scale;
  sig.dwt.Data = sig.dwt.Data * scale;
  sig.d2wt.Data = sig.d2wt.Data * scale;
  
  optimization_signals{end+1} = sig;
  if norm(sig.targ.Data) > 0
    warning(['Found nonzero target for voltage, algorithm does not ' ...
      'currently support nonzero voltage target (assumes 0).'])
  end
end

% coil currents
if gui_data.current.use_in_optimization
  sig = struct; 
  sig.name = 'ic';
  sig.calc_type = 'coil_currents';
  sig = parse_add_per_coil_waveforms(sig, gui_data.current, coils);
  sig.targ.Data = sig.targ.Data * 1e3; % kA --> A
  optimization_signals{end+1} = sig;
end

% coil current combinations
if gui_data.current_combos.use_in_optimization
  sig = struct; 
  sig.name = 'Aic';
  sig.calc_type = 'coil_current_combinations';
  sig = parse_add_weights_targs(sig, gui_data.current_combos);
  sig.targ.Data = sig.targ.Data * 1e3; % kA --> A
  sig.coil_combinations_matrix = parse_coil_combo_matrix(gui_data);
  optimization_signals{end+1} = sig;
end

% relative flux
for i = 1:3
  nam = ['rel_flux' num2str(i)];
  if gui_data.(nam).use_in_optimization
    sig = struct;
    sig.name = nam;
    sig.calc_type = 'flux_relative_multipoint';
    sig = parse_add_timeseries(sig, gui_data.(nam), {'r', 'z', 'r_multiref', 'z_multiref', 'wt_multiref'});
    sig = parse_add_weights_targs(sig, gui_data.(nam));
    optimization_signals{end+1} = sig;
  end
end

% direct flux and field 
calc_type_name_pairs = [{'flux_absolute_avg', 'flux_abs_avg'}; 
                         {'vacuum_flux_absolute_avg', 'vac_flux_abs_avg'}; 
                         {'field_absolute_radial', 'Br'}; 
                         {'field_absolute_vertical', 'Bz'}; 
                         {'vacuum_field_absolute_radial', 'Br_vac'}; 
                         {'vacuum_field_absolute_vertical', 'Bz_vac'};
                         ];
for i = 1:length(calc_type_name_pairs)
  calc_type = calc_type_name_pairs{i,1};
  nam = calc_type_name_pairs{i,2};
  if gui_data.(nam).use_in_optimization
    sig = struct;
    sig.name = nam;
    sig.calc_type = calc_type; 
    sig = parse_add_timeseries(sig, gui_data.(nam), {'r', 'z'});
    sig = parse_add_weights_targs(sig, gui_data.(nam));
    optimization_signals{end+1} = sig;
  end
end


% DEBUGGING: DO NOT MERGE
% optimization_signals{2}.wt.Data(:) = 1e-10;  % ic
% optimization_signals{3}.wt.Data(:) = 1e4;    % rel_flux1
% optimization_signals{4}.wt.Data(:) = 0;      % flux_abs_avg
% optimization_signals{4}.dwt.Data(:) = 0;
% optimization_signals{5}.wt.Data = optimization_signals{5}.wt.Data * 1e4;
% optimization_signals{6}.wt.Data = optimization_signals{6}.wt.Data * 1e4;
% gsp_inputs.settings.ic_max(:) = 400e3;
% gsp_inputs.settings.ic_min(:) = -400e3;

gsp_inputs.optimization_signals = optimization_signals;


end

% ===================
% Local functions
% ===================
function optimization_signals = remove_empty_signals(optimization_signals)
  iremove = [];
  for i = 1:length(optimization_signals)
    sig = optimization_signals{i};
    all_empty = true;
    fds = fields(sig);
    for j = 1:length(fds)
      waveform = sig.(fds{j});
      if isstruct(waveform) && isfield(waveform, 'Data')
        if ~isempty(waveform.Data) || ~isempty(waveform.Time)
          all_empty = false;
          continue
        end
      end
    end
    if all_empty
      iremove(end+1) = i;
    end
  end
  optimization_signals(iremove) = [];
end

function sig = parse_add_timeseries(sig, s, fds)
  for i = 1:length(fds)
    sig.(fds{i}).Time = s.([fds{i} '_time']);
    sig.(fds{i}).Data = s.([fds{i} '_data'])';
  end
end

function sig = parse_add_weights_targs(sig, s)
  for waveform_ = {'wt', 'dwt', 'd2wt', 'targ'}
    waveform = waveform_{:};
    sig.(waveform).Time = s.([waveform '_time']);
    sig.(waveform).Data = s.([waveform '_time_val']) * s.([waveform '_channel_val'])';
  end
end

function m = parse_coil_combo_matrix(gui_data)
  ncoils = length(fields(gui_data.ic_max));
  ncombos = numel(gui_data.current_combos.targ_channel_val);
  m = zeros(ncombos, ncoils);
  for i = 1:ncombos
    for j = 1:ncoils
      val = gui_data.current_combo_matrix.(['combo_' num2str(i-1) '_' num2str(j-1)]);
      if ~isempty(val) && ~isnan(val)
        m(i,j) = val;
      end
    end
  end
end

function fbtagcon = parse_fbtagcon_str(fbtagcon_str)
  % Parses a string like '["Ip", "Wk", "ag"]' or '['Ip', 'Wk', 'qA']' 
  % into the cell array {'Ip', 'Wk', 'ag'}  
  vars = char(strrep(sprintf(fbtagcon_str), "'", '"')); % ensure double quotes
  vars = strsplit(vars,'"'); 
  fbtagcon = vars(2:2:end); 
end

function p = interpolate_profile_to_psin(profile_struct, npsin)
  psin = linspace(0,1,npsin);
  p = struct;  
  p.Time = profile_struct.Time;  
  p.Data = interp1(profile_struct.psin, profile_struct.Data', psin)';  
end

function sig = parse_add_per_coil_waveforms(sig, s, coils)
  ncoils = length(coils);
  
  for waveform_ = {'wt', 'dwt', 'd2wt', 'targ', 'constrain'}
    waveform = waveform_{:};
    time_fds = strcat(coils, ['_' waveform '_time']);
    data_fds = strcat(coils, ['_' waveform '_data']);
  
    % each coil could be saved with different timebase, need to 
    % interpolate all onto a common timebase
    tcommon = sort(unique(struct2vec(s, time_fds)));
    
    sig.(waveform).Time = tcommon;
    sig.(waveform).Data = nan(length(tcommon), ncoils);
  
    for i = 1:ncoils
      t = s.(time_fds{i});
      data = s.(data_fds{i});
                  
      % put onto common time grid
      if ~isempty(t) && ~isempty(data) && ~all(isnan(t)) && ~all(isnan(data))               
        data = interp1(t, data, tcommon);
        
        % There is a bug in Octave where edge-value nans are not handled
        % correctly by interp1. The command "interp1([1,2,3],[1,2,nan], 2)"
        % incorrectly produces the result nan. This is hack to paint in
        % nans, so that they will not result in nans for future interp1 
        % operations. Method = 4 corresponds to using a nearest neighbor 
        % result for the edge-values (inpaint_nans([1,2,nan], 4) produces [1,2,2])         
        method = 4;
        sig.(waveform).Data(:,i) = inpaint_nans(data, method);

      end
    end
  end
end
