% USAGE:    
%   plasma_params = define_plasma_params(settings)
%   
% DESCRIPTION: 
%  Define core plasma properties of the equilibrium evolution. 
%
% INPUTS: 
%   settings - settings struct
%
% OUTPUTS: 
%   plasma_params - struct with waveforms corresponding to the evolution of
%                   the core plasma properties. 
%
% ADDITIONAL INFO: 
%
% Several core plasma properties must be specified in order to constrain
% the Grad-Shafranov solve P' and FF' profiles, and also the flux evolution 
% when using the Ejima model. The required plasma_params fields are: 
%
% (Rp.Time, Rp.Data) - plasma resistance vs time, only required if
%    settings.specify_psibry_mode = 'ejima'
% 
% (pprime.Time, pprime.Data) - pprime profile shape vs time, on a uniform
% psiN basis. 
%
% [WHEN USING SETTINGS.PICARD_ALGO = 'GSPULSE']
% 
% In this case, using only 1 ttprime basis function is supported, which
% should be specified in (ttprime.Time, ttprime.Data)
%
% The 'gspulse' method works by constraining plasma current (Ip) and total
% thermal energy (Wk). In this case, specify these via: (Ip.Time, Ip.Data)
% and (Wk.Time, Wk.Data)
%
% [WHEN USING SETTINGS.PICARD_ALGO = 'FBT']
%
% FBT supports many options for constraining the plasma picard iteration. 
% To see how GSPulse interacts with FBT, see fbt_update_psipla.m. These
% depend on the exact solve type and settings in L. GSPulse supports a
% subset of these which is L.P.bfct = @bfabmex and L.P.bfct = @bf3imex. 
% 
% When using L.P.bfct = @bfabmex, the basis functions do not need to be
% specified (calculated internally), and the user just needs to specify
% plasma current (Ip.Time, Ip.Data) and thermal energy (Wk.Time, Wk.Data). 
%
% When using L.P.bfct = @bf3imex, specify 3 basis functions vs time using
% (pprime.Time, pprime.Data), (ttprime1.Time, ttprime1.Data), and
% (ttprime2.Time, ttprime2.Data). 
%
% There are multiple constraints that can be implemented. A typical example
% is the plasma current, thermal energy, and q on axis (L.P.fbtagcon =
% {'Ip', 'Wk', 'qA'}. In this case, the user must specify values for each
% of these constraints (Ip.Time, Ip.Data), (Wk.Time, Wk.Data), and
% (qA.Time, qA.Data). Additionally, the qA constraint needs to specify the
% plasma toroidal field and radius (rBt.Time, rBt.Data). 
%
% In general, there are many options for constraining FBT which are set by 
% flags in the L object, and then whatever parameters are needed for the 
% FBT LX object must be specified as waveforms here in plasma_params. 
