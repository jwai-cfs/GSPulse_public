% USAGE:    
%   settings = define_init()
%
% DESCRIPTION: 
%
% Define the initial condition for the optimization. These will be enforced
% as equality constraints. Use nan to indicate no constraint.  
%
% INPUTS: 
%   None
%
% OUTPUTS: 
%  init - struct with fields 'x1', 'u0', and 'u1'
%
%  init.x1 - vector of coil currents and vessel currents, corresponding to 
%            t=1 (the first equilibrium that is solved for). If using 
%            eigenmode approach to compress the vessels elements, then the 
%            currents should be provided as the eigenmode currents. Should 
%            have dimension (# of coils + # of vessel eigenmodes x 1), or 
%            if no constraints set init.x1=nan. Individual elements of the 
%            x1 vector can be set to nan to specify no constraint on that
%            element. 
%        
%  init.u0 - vector of coil voltages corresponding to t=0, (the voltages
%            at the time step prior to the first equilibrium). These 
%            voltages do not affect the dynamic evolution of the
%            trajectory, however, they can influence the optimization if
%            there are time-derivative or time-second-derivative weights on
%            the voltage. Should either have dimension (# coils x 1) or 
%            set init.u0 = nan. 
% 
%  init.u1 - vector of coil voltages corresponding to t=1. 
