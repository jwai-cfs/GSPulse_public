% Define various settings for the optimization. Try 'help
% _define_settings.m' for details. 

function settings = define_settings()

s = struct; 
s.plotlevel = 1; 
s.interval_t{1} = 0.05:0.02:1; 
s.niter = 5;  
s.verbose = 1;
s.picard_algo = 'gspulse';
s.do_final_boundary_trace = 0;
s.use_spline_basis = 1;
s.spline_basis_ratio = 0.7;
s.nvessmodes = 40;      
s.specify_psibry_mode = 'ejima'; 
s.calc_strike_pts = 0;
s.dosave = 0;
s.qpsolver = 'quadprog';

s.fds2control = {};
s.fds2control{end+1} = 'ic';                   % current in the PF coils
s.fds2control{end+1} = 'diff_psicp_psix1';     % flux difference between the shape control points and x-point 1
s.fds2control{end+1} = 'diff_psicp_psitouch';  % flux difference between the shape control points and the touch point
s.fds2control{end+1} = 'psibry';               % value of the boundary flux
s.fds2control{end+1} = 'psix_r';               % gradient of the flux w.r.t. r measured at the x-points
s.fds2control{end+1} = 'psix_z';               % gradient of the flux w.r.t. z measured at the x-points

s.inject_model_offset = 0;

s.enforce_voltage_limits = 1;
s.vmax = [4048 1012 2024 2024 3036 2024 2024 1012]';
s.vmin = -s.vmax;

s.enforce_current_limits = 1;
s.ic_max = [20  15 15 8 0 8 15 15]' * 1e3;
s.ic_min = [-20 0 0 -13 -24 -13 0 0]' * 1e3;

settings = s;