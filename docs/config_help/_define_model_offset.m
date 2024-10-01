% USAGE:    
%   model_offset = define_model_offset(tok)
%
% DESCRIPTION: 
%   Prescribe an estimate of the modelling errors which will GSPulse will
%   attempt to compensate for
%
% INPUTS: 
%  tok - the tokamak geometry object
%
% OUTPUTS: 
%  model_offset - struct with waveform for 'psizr'
%  (model_offset.psizr.Time, model_offset.psizr.Data) - waveform describing 
%     the estimated modelling error of the flux on grid vs time
%
% ADDITIONAL INFO
% 
% For theoretical studies, the model offset should generally be 0. The use
% case here is for the scenario where:
%
%  - a GSPulse scenario was created
%  - an experiment was run, revealing differences between the actual
%    tokamak and the GSPulse predictions
%  - now run GSPulse again, and attempt to compensate for the modelling
%    discrepancies
%
% As an example, let us consider the error of achieving some parameter.
% The total error is:
%
% e_tot = e_modelled + e_unmmodelled
%
% That is the total error is the sum of modelled error and un-modelled
% errors. GSPulse attempts to minimize the modelled error using its own 
% internal models. However, if we have some estimate of the e_ummodelled 
% then we can inject this into the GSPulse optimization and minimize e_tot
% instead. 
%
% The unmodelled error should be specified as the total flux on the grid
% vs. time, from which specific flux error differences and fields can be
% interpolated. 
