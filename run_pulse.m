function [soln, gsp_inputs] = run_pulse(tokamak, pulseid, inputs_only, doplot, varargin)
% USAGE:    
%   [soln, gsp_inputs] = run_pulse(tokamak, pulseid, inputs_only, doplot,...
%                            'setting1',value,'setting2',value)
%
% DESCRIPTION: 
%
% Run GSPulse for a specific set of config files
%
% INPUTS: 
%  tokamak - string, name of the tokamak for which to run the pulse
%  pulseid - int, the config files for this pulse are assumed
%            to be located in 'tokamaks/<tokamak>/pulses/<pulseid>'
%  inputs_only - bool (default false), if false, runs the config files to 
%                generate the gsp_inputs and then runs GSPulse. If false, 
%                only generates the gsp_inputs. 
%  doplot - bool (default true), whether to plot the GSPulse results
%  [optional]: setting,value pairs to override settings from cmd line
%
% OUTPUTS: 
%  soln - GSPulse solution struct with fields
%    t: timebase of solve
%    eqs: equilibrium description at each time
%    signals: waveforms for all of the signals that the optimizer computes
%             such as the coil currents, voltages, and flux errors
%  gsp_inputs - struct containing all of the inputs needed to execute a 
%               GSPulse run
%
% ADDITIONAL INFO
%
% Each run of GSPulse is created by creating a set of config files and
% placing these in the 'tokamaks/<tokamak>/pulses/<pulseid>'. To create a
% new pulse, the easiest way is to copy an existing pulse and modify the
% content as needed. 
%
% The required config files for a GSPulse run are:
% * define_settings.m
% * define_tok.m
% * define_init.m
% * define_shapes.m
% * define_plasma_params.m
% * define_targets.m
% * define_weights.m
%
% Help for each of these files can be found in docs/config_help or by
% running 'help _<config file>' (note the underscore): For example,
%
% help _define_general_settings.m
%
% A written description of the GSPulse algorithm is available in
% docs/gspulse_algorithm.pdf

if ~exist('inputs_only', 'var'), inputs_only = 0; end
if ~exist('doplot', 'var'), doplot = 1; end

user_dir = pwd();
cd([getenv('GSROOT') '/tokamaks/' lower(tokamak) '/pulses/' num2str(pulseid)]);

settings = define_settings();

% overwrite user-defined parameters
for ii=1:2:numel(varargin)
  fieldname = varargin{ii};
  settings.(fieldname) = varargin{ii+1};
end

% check that settings has enough fields specified
settings = check_default_settings(settings, tokamak);                               

switch settings.picard_algo
  case 'fbt'
    [tok, L] = define_tok(settings);    
  case 'gspulse'    
    tok = define_tok(settings);
    L = struct();
end

init             = define_init();
shapes           = define_shapes(settings, tok);
plasma_params    = define_plasma_params(settings);
targs            = define_targets(settings, shapes, tok);
weights          = define_weights(settings, targs, plasma_params, tok);

if settings.inject_model_offset
  model_offset = define_model_offset(tok);
else
  model_offset = struct;
end

gsp_inputs = variables2struct(settings, L, tok, shapes, plasma_params, ...
  targs, weights, init, model_offset);

if settings.dosave
  fname = sprintf('gsp_inputs_%d.mat',pulseid);
  fprintf('saving %s',fname);
  gg = rmfield(gsp_inputs,'L');
  gg.G = gsp_inputs.L.G;
  gg.P = gsp_inputs.L.P;
  save(fname,'-struct','gg');
end

if ~inputs_only
  soln = GSPulse(gsp_inputs);
  if doplot 
    plot_soln(soln, gsp_inputs);
  end
else
  soln = [];
end

cd(user_dir)
