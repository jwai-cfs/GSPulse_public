% Define time-dependent weights for each of the settings.fds2control (and
% also the actuator voltage). 
%
% weights.wts is the weight on the value of the parameter. weights.dwts is
% the weight on the derivative of the value, and d2wts is the weights on the 
% second derivative. Parameters are not normalized so the weights may need
% to span orders of magnitude depending on the normal range of each 
% measured variable. 
%  
% To view a summary plot of the weights, set settings.plotlevel >= 2

function weights = define_optimization_weights(settings, targs, ~, tok)

% ==================
% Initialize weights
% ==================

% timebase on which to define all weights
nt = 101;
t = linspace(0, 1, nt)';   

% initialize weights 
[wts, dwts, d2wts] = initialize_weights_to_zero(t, settings, targs, tok);

% define some dimensions 
ncp = size(wts.diff_psicp_psix1.Data, 2);  % number of shape control points

% ============================
% Custom-populate the weights
% ============================

% weight on psibry (for surface voltage and driving Ip)
wts.psibry.Data(:) = 5e6;  

% wt on flux err vs touch-point starts on and turns off as plasma diverts
wts.diff_psicp_psitouch.Data = sigmoidn(t, 0.15, 0.2, 1, 0) * ones(1,ncp) * 5e6;

% wt on flux err vs x-point starts off and turns on as plasma diverts
wts.diff_psicp_psix.Data  = sigmoidn(t, 0.15, 0.2, 0, 1) * ones(1,ncp) * 5e6;

% weight on flux gradient turns on as plasma diverts
wts.psix_r.Data(:) = sigmoidn(t, 0.15, 0.25, 0, 1) * 3e7;
wts.psix_z.Data(:) = sigmoidn(t, 0.15, 0.25, 0, 1) * 3e7;

% weight the outer boundary point even higher
wts.diff_psicp_psitouch.Data(:) = wts.diff_psicp_psitouch.Data * 30;
wts.diff_psicp_psix1.Data(:)     = wts.diff_psicp_psix.Data * 30;

% no weight on absolute voltage
wts.v.Data = zeros(nt,tok.nc);   

% weight on the change in voltage - as is, the weights are roughly
% proportional to 1/coil inductance, but not too sensitive to the weighting
dwts.v.Data = ones(nt,1) * [1 1 1 1 0.2 1 1 1];


% ==============
% Save and plot 
% ==============
weights = variables2struct(wts, dwts, d2wts);

% plot weights
if settings.plotlevel >= 2
  fds = [settings.fds2control(:); 'v'];

  h = plot_structts(wts, fds, [], 3, [], 'linewidth', 1.5);
  sgtitle('weights.wts', 'fontsize', 14); drawnow

  h = plot_structts(dwts, fds, [], 3, [], 'linewidth', 1.5);
  sgtitle('weights.dwts', 'fontsize', 14); drawnow

  h = plot_structts(d2wts, fds, [], 3, [], 'linewidth', 1.5);
  sgtitle('weights.d2wts', 'fontsize', 14); drawnow
end








