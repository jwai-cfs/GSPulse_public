function [soln, gsp_inputs] = run_pulse(tokamak, pulse_id, varargin)
% =========================================================================
% EXAMPLE (basic):
% tokamak = 'sparc';
% pulse_id = 113;
% [soln, gsp_inputs] = run_pulse(tokamak, pulse_id);
%
% EXAMPLE (switch qpsolver from default): 
% settings_overrides = struct('qpsolver', 'scs');
% [soln, gsp_inputs] = run_pulse(tokamak, pulse_id, ...
%     'settings_overrides', settings_overrides);
% =========================================================================

% parse arguments (Note: Octave does not support the 'arguments' syntax, 
% use inputParser syntax, supported by both Matlab & Octave)
p = inputParser;
addParameter(p, 'gui_data_fp', 'default');
addParameter(p, 'inputs_only', 0);
addParameter(p, 'gsp_inputs', []);
addParameter(p, 'settings_overrides', struct);
parse(p, varargin{:});
options = p.Results;

% define the GUI data filepath
if isnumeric(pulse_id)
  id_str = num2str(pulse_id);
else
  id_str = pulse_id;
end
if strcmp(options.gui_data_fp, 'default')  
  options.gui_data_fp = [getenv('GSROOT'), '/tokamaks/', tokamak, '/runs/processed_inputs/', id_str, '_gui_state_processed.json'];
end

% parse gui data and form gspulse inputs
if isempty(options.gsp_inputs)
  gsp_inputs = parse_gsp_inputs_from_gui(options.gui_data_fp, options.settings_overrides);
  gsp_inputs = validate_gspulse_input_config(gsp_inputs);
else
  gsp_inputs = options.gsp_inputs;
end

% run GSPulse
if ~options.inputs_only
  soln = GSPulse(gsp_inputs);
  if gsp_inputs.settings.plotlevel>0 % visualize solution
    plot_soln(soln, gsp_inputs);
  end
else
  soln = [];
end
  