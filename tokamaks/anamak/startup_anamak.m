function startup_anamak()
% Do any device-specific startup tasks, such as adding device paths to the 
% Matlab path. 

GSROOT = getenv('GSROOT');
addpath(fullfile(GSROOT, 'tokamaks/anamak/device'));