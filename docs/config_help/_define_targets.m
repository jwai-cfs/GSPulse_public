% USAGE:    
%   targs = define_targets(settings, shapes, tok)
%
% DESCRIPTION: 
%   Define targets for the GSPulse optimization
%
% INPUTS: 
%  settings - settings struct, see help _define_settings.m
%  shapes   - shapes   struct, see help _define_shapes.m
%  tok      - tok      struct, see help _define_tok.m
%
% OUTPUTS: 
%   targs - targs struct, with waveforms defined for each of the fields in 
%           settings.fds2control, as well as the voltage (targs.v)
%
% ADDITIONAL INFO:
%  
%   settings.fds2control specifies fields will be explicitly optimized for.
%   For example, if fds2control{1} = 'ic' (that is, the coil currents),
%   then the optimizer will attempt to have the equilibrium coil currents
%   match the waveform of (targs.ic.Time, targs.ic.Data) and penalize any
%   deviation with weights corresponding to (weights.wts.ic,
%   weights.dwts.ic, weights.d2wts.ic) 
%
%   It is often the case that the target is zero, for example, we often
%   want to control flux differences between the target boundary and the
%   x-point to be zero. However, this is not always the case, such as if it
%   is desired to specify that a point lies slightly off of the boundary,
%   or specifying a drsep value. 
