function [soln, gsp_inputs] = plot_saved_run(tokamak, pulse_id, options)
% =========================================================================
% EXAMPLE:
% tokamak = 'sparc';
% pulse_id = 113;
% [soln, gsp_inputs] = plot_saved_run(tokamak, pulse_id);
% =========================================================================
arguments
  tokamak
  pulse_id
  options.gui_data_fp = 'default';
  options.soln_fp = 'default';
end

% parse inputs
if isnumeric(pulse_id)
  id_str = num2str(pulse_id);
else
  id_str = pulse_id;
end
if strcmp(options.gui_data_fp, 'default')  
  options.gui_data_fp = [getenv('GSROOT'), '/tokamaks/', tokamak, '/runs/processed_inputs/', id_str, '_gui_state_processed.json'];
end
if strcmp(options.soln_fp, 'default')
  options.soln_fp = [getenv('GSROOT'), '/tokamaks/', tokamak, '/runs/soln/', id_str, '_soln.mat'];
end

% parse gui data and form gspulse inputs
gsp_inputs = parse_gsp_inputs_from_gui(options.gui_data_fp);
gsp_inputs = validate_gspulse_input_config(gsp_inputs);

% load solution and plot
soln = load(options.soln_fp).soln;
plot_soln(soln, gsp_inputs);
  