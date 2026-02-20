function startup_gspulse(tokamak, varargin)
% Parse arguments. (Note: Octave does not support the 'arguments' syntax, 
% use inputParser syntax, supported by both Matlab & Octave)
p = inputParser;
addParameter(p, 'meqpath', 'default');
parse(p, varargin{:});
options = p.Results;

% Set GSPulse root
current_dir = pwd();
GSROOT = fileparts(fileparts(fileparts(fullfile(mfilename('fullpath'))))); 
setenv('GSROOT', GSROOT);
addpath(genpath(fullfile(GSROOT, 'src/matlab')))

% Set MEQ location
is_octave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if strcmp(options.meqpath, 'default')
  if is_octave
    options.meqpath = fullfile(GSROOT, 'build', 'meq_octave_build');
  else % matlab
    options.meqpath = fullfile(GSROOT, 'build', 'meq_matlab_build');
  end
end

% Add MEQ paths
if ~isfolder(options.meqpath)
  warning('Specified meqpath does not exist: %s', options.meqpath);
else
  fprintf('Using meqpath @ %s\n', options.meqpath)
  if is_octave
    warning('off', 'Octave:shadowed-function')
    addpath(fullfile(options.meqpath, '_octave'))
    addpath(fullfile(GSROOT, 'src/octave'));
  end
  addpath(options.meqpath)
  addpath(fullfile(options.meqpath, 'IDS'))

  % If not added here, these get added the first time 'fbt' is called, but
  % preference is for all paths to be added at startup time for visibility.
  addpath(fullfile(options.meqpath, 'genlib'))  
  addpath(fullfile(options.meqpath, 'solver'))
end

% Check that MEQ is compiled
meq_is_compiled = exist('greenemmex', 'file') == 3;
if ~meq_is_compiled
  error(['Detected that MEQ @ %s has not been compiled.' ...
    'Have you run the install script at "scripts/install.sh"? ' ...
    'See README.md for help.'], options.meqpath)
end

% Clear previous devicepaths
warnstate = warning('off', 'all');  % suppress warnings for the next command
rmpath(genpath(fullfile(GSROOT, 'tokamaks')))
rmpath(genpath(fullfile(GSROOT, 'submodules/MEQ-CFS/device')))
warning(warnstate);  % warnings back on


% Perform device-specific startup tasks
% assumes there is a file called "add_<tokamak>_startup_paths.m"
mfile = sprintf('%s/tokamaks/%s/startup_%s.m', GSROOT, tokamak, tokamak);
if isfile(mfile)
  run(mfile);
else
  warning('Could not find device startup file:\n%s', mfile);
end

% Compile SCS for Matlab
if ~is_octave 
  disp('Compiling SCS ...')
  addpath(fullfile(GSROOT, 'submodules/scs-matlab'));
  if ~(exist('scs_direct', 'file') == 3) 
    cd(fullfile(GSROOT, 'submodules/scs-matlab'))
    make_scs();
  end
end

cd(current_dir)
disp('All GSPulse startup complete.')
