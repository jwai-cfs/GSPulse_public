% USAGE:    
%   settings = define_general_settings()
%
% DESCRIPTION: 
%   Define top-level settings for a GSPulse run. 
%
% INPUTS: 
%   None
%
% OUTPUTS: 
%   settings - struct, with values set for the following required or
%              optional fields
%
% ADDITIONAL INFO:
%
% Each pulse must have a version of this script to define top-level solver
% parameters such as the time-basis on which to solve and the number of 
% iterations to perform. The required and optional flags that must be
% defined are: 
%
% [REQUIRED PARAMETERS]
% 
% interval_t: 
%
% cell array, with each each element of the cell containing a
% vector of solution times. This defines the timebase(s) for solving the
% optimization and equilibria. As an example: 
%
%   interval_t{1} = 0.2: 0.2: 5;
%   interval_t{2} = 4.9: 0.1: 7;
%   interval_t{3} = 6  : 0.1: 10;
%
% would solve for the solution in 3 different stages (e.g. the first 
% stage solves for equilibria starting from 0.2 seconds to 5 seconds, on 
% a 0.2s interval). The time basis for each stage can have different
% spacing, but each interval should have uniform spacing. In order to
% have smooth transitions across stages, it is necessary to overlap the
% endpoints of each stage slightly (recommended, 2 samples overlap). In 
% this example, the second stage is overlapping the first stage by 2
% samples (i.e. at 4.9 and 5.0 sec). The third stage is overlapping the 
% second stage by many more samples than strictly necessary which has
% higher computational cost, but may also give better trajectories.
%
% fds2control: 
%
% cellarray, with each element containing a string corresponding to the
% variables that will be explicitly controlled by the optimization
% algorithm. For every variable that is included, the user must set 
% corresponding weights and targets, and there must be a method for 
% measuring the variable (see measure_ys.m) and defining the linearization
% (see output_model.m). The currently-supported options are: 
%
%  ic  - current in the PF coils
%  Aic - linear combinations of PF coil currents, the linear combinations
%        are speficied by defining tok.Aic
%  psix_r - flux gradient wrt R at the target xpts (shapes.rx, shapes.zx)
%  psix_z - flux gradient wrt Z at the target xpts (shapes.rx, shapes.zx)
%  psibry - flux value at the target boundary (shapes.rbdef, shapes.zbdef)
%  diff_psicp_psix1 - the difference in flux between the target control
%    points (shapes.rb, shapes.zb) and the first x-point
%  diff_psicp_psix2 - the difference in flux between the target control
%    points (shapes.rb, shapes.zb) and the second x-point
%  diff_psicp_psitouch - the difference in flux between the target control
%    points (shapes.rb, shapes.zb) and the target touch point
%    (shapes.rtouch, shapes.ztouch)
%  diff_psix1_psix2 - difference in flux between the first and second
%    xpts (e.g. for controlling drsep)
%  diff_psix2_psix3 - difference in flux between the second and third
%    xpts (e.g. for controlling drsep). One might have 3 xpoints for
%    example in an X-Point Target divertor. 
%  diff_psisp_psix1 - difference in flux between the target strike points
%    (shapes.rstrike, shapes.zstrike) and the first x-point.
%  diff_psisp_psix2 - difference in flux between the target strike points
%    (shapes.rstrike, shapes.zstrike) and the second x-point.
%
% vmax: 
% maximum voltage for each coil, must have dimensions [# coils x 1]
%
% vmin: 
% minimum voltage for each coil, must have dimensions [# coils x 1]
%
% ic_max:
% maximum current allowed in each coil, must have dimensions [# coils x 1]
% 
% ic_min
% minimum current allowed in each coil, must have dimensions [# coils x 1]
% 
% 
% [OPTIONAL PARAMETERS]
% 
% plotlevel (default, 1):
% int, with values 0, 1, or 2. 0 = don't make any plots, 1 = plot the 
% output results, 2 = plot inputs and settings and outputs
% 
% verbose (default, 1):
% int, with values 0 (minimal) or 1 (verbose)
% 
% niter (default, 10):
% int, number of Grad-Shafranov Picard iterations to perform 
% 
% use_spline_basis (default, False):
% bool, whether to perform a spline basis compression when 
% solving the MPC quadratic program
% 
% spline_basis_ratio (default, 0.7):
% float, if using a spline basis, a target compression ratio for how much to 
% compress the optimization. E.g., if there are originally N points in the
% interval, the QP will solve for ~spline_basis_ratio*N spline knot points
% 
% nvessmodes (default, 40):
% number of vessel modes to include. Must be consistent with the dimensions 
% of the resistances and mutuals in the "tok" object
%
% specify_psibry_mode (default, 'direct'): 
% Either 'ejima' or 'direct'. If 'direct', then the target boundary flux
% (psibry) is specified as a user-waveform in targs.psibry. If 'ejima',
% then the user must specify a target waveform for the plasma resistance in
% plasma_params.Rp, and then GSPulse will compute the corresponding
% targs.psibry as a function of the resistance and internal inductance (see
% compute_target_psibry.m). 
% 
% enforce_voltage_limits (default, True):
% bool, whether or not to enforce the voltage limits (specified in 
% settings.vmax, settings.vmin) 
%
% enforce_current_limits (default, False):
% bool, whether or not to enforce the voltage limits (specified in 
% settings.ic_max, settings.ic_min) 
%
% picard_algo (default, 'fbt'):
% Either 'fbt' or 'gspulse'. Using 'fbt' requires access to the MEQ
% submodule but is a little bit more stable. This method calls out to the
% picard iteration in the FBT code part of the MEQ suite. It requires the 
% input geometry to be specified in the MEQ format (called "L" in FBT), and
% then the equilibria will include all of the relevant FBT output formatting 
% (called "LY"). 
% 
% If the option is 'gspulse', then GSPulse will use the built-in picard
% iteration solver and boundary-finding methods. With this option, the
% input geometry is specified via the 'tok' object (not the MEQ 'L'
% object). See also tok.m, L2tok.m, tok_data_struct2tok.m
%
% do_final_boundary_trace (default, False):
% bool, whether or not to trace the plasma boundary on the last iteration.
% The traced boundary is returned in [soln.eqs{}.rbbbs, soln.eqs{}.zbbbs]
%
% inject_model_offset (default, False):
% bool, whether or not to inject an approximation of the modelling error
% into the solver. See help _define_model_offset.m. 
