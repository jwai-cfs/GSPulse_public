% USAGE:  
%   weights = define_weights(settings, targs, plasma_params, tok)
%
% DESCRIPTION: 
%   Define weights for the GSPulse optimization
%
% INPUTS: 
%  settings - settings struct, see help _define_settings.m
%  targs    - targs    struct, see help _define_targs.m
%  shapes   - shapes   struct, see help _define_shapes.m
%  tok      - tok      struct, see help _define_tok.m
%
% OUTPUTS: 
%  weights - struct with fields, 'wts', 'dwts', 'd2wts'
%  weights.wts  - struct with waveforms weights corresponding to each of
%                 the settings.fds2control as well the power supply
%                 voltages ('v'). 'wts' is the weight of the value of that
%                 parameter. 
%  weights.dwts - sames as .wts, but corresponding to the derivative w.r.t
%                 time for that parameter. Used for penalizing ramping
%                 trajectories.
%
%  weights.d2wts - sames as .wts, but corresponding to the 2nd derivative 
%                  w.r.t time for that parameter. Used for penalizing
%                  non-smooth trajectories. 
%
% ADDITIONAL INFO:
%  
% Example: if settings.fds2control = {'ic', 'diff_psicp_psix1'}, that is,
% coil currents and flux differences between the target boundary and
% x-point, then weights will contain:
%
% (weights.wts.ic.Time, weights.wts.ic.Data) - waveform for the currents
% (weights.wts.diff_psicp_psix1.Time, weights.wts.diff_psicp_psix1.Data) -
%   waveform for the flux errors
% (weights.wts.v.Time, weights.wts.v.Data) - waveform for the voltage 
% 
% and corresponding waveforms for the derivatives (dwts) and
% second-derivatives (d2wts). 
% 
% Note that the various parameters are NOT normalized, so the weights may 
% need to span orders of magnitude depending on the normal range of each 
% measured variable. 

