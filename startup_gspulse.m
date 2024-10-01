% set environment variable
GSROOT = fileparts(mfilename('fullpath'));
setenv('GSROOT', GSROOT);

% add path
addpath(GSROOT)
addpath(genpath([GSROOT '/src']))
addpath([GSROOT '/MEQ'])
addpath([GSROOT '/MEQ/IDS'])
addpath([GSROOT '/docs/config_help'])

isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if isOctave
  fprintf('adding MEQ/_octave path for Octave compatibility\n')
  addpath([GSROOT,'/MEQ/_octave'])
end

fprintf("startup_gspulse.m: Added paths. For help, see README.md or try 'help run_pulse.m'\n")
