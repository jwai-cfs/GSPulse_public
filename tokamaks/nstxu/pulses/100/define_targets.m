% Define time-dependent targets for each of the settings.fds2control 
%
% Most of these are zero. In each case the error to be minimized is this
% target value minus the value measured at the corresponding equilibrium. 
%
% To view a summary plot of the targets, set settings.plotlevel>=2

function targs = define_optimization_targets(settings, shapes, tok)

% define dimensions
nt = 101;
t = linspace(0, 1, nt)';
ncp = size(shapes.rb.Data,2); % number of shape control pts
nxp = 1;                      % number of x-points
nsp = 2;                      % number of strike points
nc = tok.nc;                  % number of coils

% define values for isoflux targets
% .................................

% flux gradient at target x-point(s)
targs.psix_r.Data = zeros(nt,nxp);
targs.psix_r.Time = t;

% flux gradient at target x-point(s)
targs.psix_z.Data = zeros(nt,nxp);
targs.psix_z.Time = t;

% flux error - control points vs x-point
targs.diff_psicp_psix1.Data = zeros(nt,ncp);
targs.diff_psicp_psix1.Time = t;

targs.diff_psicp_psix2.Data = zeros(nt,ncp);
targs.diff_psicp_psix2.Time = t;

% flux error - strike points vs x-point
targs.diff_psisp_psix1.Data = zeros(nt,nsp);
targs.diff_psisp_psix1.Time = t;

targs.diff_psisp_psix2.Data = zeros(nt,nsp);
targs.diff_psisp_psix2.Time = t;

% flux error - control points vs touch point
targs.diff_psicp_psitouch.Data = zeros(nt,ncp);
targs.diff_psicp_psitouch.Time = t;

% drsep (in flux-space)
targs.diff_psix1_psix2.Data = zeros(nt,1);
targs.diff_psix1_psix2.Time = t;

% coil currents
targs.ic.Data = zeros(nt,nc);
targs.ic.Time = t;

% if settings.specify_psibry_directly == 0, psibry target will be computed 
% automatically to satisfy the Ejima model via compute_target_psibry.m. 
% Otherwise, the target psibry should be specified here. 
targs.psibry.Time = t;
targs.psibry.Data = t*nan;

% transpose data if necessary
targs = check_structts_dims(targs);

% plot timetraces of the targets
if settings.plotlevel >= 2
  h = plot_structts(targs, settings.fds2control, [], 3, [], 'linewidth', 1.5);
  sgtitle('targs', 'fontsize', 14); drawnow
end









