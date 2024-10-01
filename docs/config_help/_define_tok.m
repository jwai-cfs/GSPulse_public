% USAGE:    
%   tok = define_tok(settings);      % if settings.picard_algo = 'gspulse'
%   [tok, L] = define_tok(settings)  % if settings.picard_algo = 'fbt'
%   
% DESCRIPTION: 
%  Define the geometry description of the tokamak. 
%
% INPUTS: 
%   settings - settings struct, see help _define_settings.m
%
% OUTPUTS: 
%   tok - tokamak geometry object with mutual inductance and resistance
%         tables. See description below. Note that GSPulse supports
%         creating this object from the geometry objects of TokSys and/or
%         MEQ. For conversion from TokSys, see tok_data_struct2tok.m. For
%         conversion from MEQ, see L2tok.m
%   L   - the MEQ geometry and solver object, used when calling the FBT 
%         algorithm Picard solver. 
%
% ADDITIONAL INFO: 
%
% tok contains the geometry information, mutual inductance tables, and 
% coil resistances needed for defining the dynamic system and solving the
% GS equation. 
% 
% In general, the mutuals and resistances should represent the
% *circuit-connected* mutuals and resistances. (e.g., if a PF coil consists 
% of two elements connected in series, then the mutual inductance should
% reprsent the contribution from both those elements). Circuit connections
% are not specified elsewhere, if wanting to evaluate different circuit
% connections then multiple tok objects must be created. 
%
% GSPulse supports creating this object from the geometry objects of 
% TokSys and/or MEQ. For conversion from TokSys, see tok_data_struct2tok.m.
% For conversion from MEQ, see L2tok.m. 
% 
% The tok object must contain the following fields: 
%
% rg - radial grid
% zg - vertical grid
% nr - number of radial grid points
% nz - number of vertical grid points
% nc - number of coil elements (if using circuits, should correspond to 
%      the number of coil circuits)
% nv - number of vessel elements (if using vessel modes, should correspond
%      to the number of vessel modes)
% rl - R of the limiter
% zl - Z of the limiter
% resc - coil resistances
% resv - vessel resistances
% mcc - coil-to-coil mutuals, dimension [nc x nc]
% mcv - coil-to-vessel mutuals, dimension [nc x nv]
% mvv - vessel-to-vessel mutuals, dimension [nv x nv]
% mpc - plasma grid to coil mutuals, dimension [nz*nr x nc]
% mpv - plasma grid to vessel mutuals, dimension [nz*nr x nv]
% mpp - plasma grid to plasma grid mutuals, dimension [nz*nr x nz*nr]
% ccnames - cell array of coil names
% (rgg, zgg) - meshgrid of (rg,zg)
% outside_vessel_mask - flag of which points lie outside the limiter,
%    dimension [nz x nr]