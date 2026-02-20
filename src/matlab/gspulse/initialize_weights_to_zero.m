function [wts, dwts, d2wts] = initialize_weights_to_zero(t, settings, targs, tok)
% =========================================================================
% Description: 
%  initializes all of the absolute weights (wts), weights on the
%  derivatives (dwts), and weights on the second derivatives (d2wts) to
%  zeros of the appropriate data dimension. 
%
% Inputs: 
%  t - timebase used for initializing weight waveforms
%  settings      - settings struct, see help _define_settings.m
%  targs         - targets struct, see help _define_targets.m
%  tok           - tokamak geometry struct, see help _define_tok.m
%
% Outputs: 
%  wts, dwts, d2wts - struct of timeseries objects, with weight values
%                     initialized to zero
%
% =========================================================================
cv = define_data_indices(settings, targs, tok);
N = length(t); 

% initialize target weights
for dum = cv.ynames(:)'
  fd = dum{:};
  wts.(fd).Time = t(:);
  dwts.(fd).Time = t(:);
  d2wts.(fd).Time = t(:);

  ny = size(targs.(fd).Data, 2);
  wts.(fd).Data  = zeros(N,ny);    
  dwts.(fd).Data = zeros(N,ny); 
  d2wts.(fd).Data = zeros(N,ny);

end

% initialize voltage weights
wts.v.Time = t(:);
dwts.v.Time = t(:);
d2wts.v.Time = t(:);
wts.v.Data = zeros(N,tok.nc);
dwts.v.Data = zeros(N,tok.nc);
d2wts.v.Data = zeros(N,tok.nc);