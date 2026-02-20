function [soln, gsp_inputs] = plot_python_soln(tokamak, pulse_id)
% Wrapper function for importing data saved from python and plotting

config_dir = fullfile(getenv('GSROOT'), 'tokamaks', tokamak, 'pulses', num2str(pulse_id));
soln_fp = fullfile(config_dir, 'soln.mat');
gspulse_inputs_pyconfig_fp = fullfile(config_dir, 'gspulse_inputs_pyconfig.mat');

fprintf('Loading soln from %s\n', soln_fp);
soln = load(soln_fp).soln;

fprintf('Loading inputs from %s\n', gspulse_inputs_pyconfig_fp);
gspulse_inputs_pyconfig = load(gspulse_inputs_pyconfig_fp).gspulse_inputs_pyconfig;

gsp_inputs = gsp_inputs_py2mat(gspulse_inputs_pyconfig);
gsp_inputs = validate_gspulse_input_config(gsp_inputs);

plot_soln(soln, gsp_inputs)
