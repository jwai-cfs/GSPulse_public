function init = define_init_from_prev_stage(timedat, prevstagesoln, cv, tok)
% =========================================================================
% Description: 
%  when using multiple intervals to design/optimize the trajectory, the
%  initial condition for the next stage must overlap with final condition
%  from the previous stage, to ensure continuity. This function extracts the
%  appropriate initial conditions that are needed for the next stage from
%  the previous stage. 
%
% 
% Inputs: 
%  timedat - struct with timing information about the next stage
%  prevstagesoln - struct with the GSPulse solution output from the
%                   previous stage
%  cv - struct with information on the databus indices and grouping, see 
%       define_data_indices.m
%  tok - tokamak geometry object
% 
% Outputs: 
%  init - object with initial condition data for the next stage
%
% Additional info: 
%   Appendix B.2 of the gspulse_algorithm.pdf describes optimization across
%   multiple intervals in more detail. 
%
% =========================================================================
s = prevstagesoln;
s.mpcsoln.pcurrt = s.pcurrt;

dt = timedat.dt;
t0 = timedat.t0;
t1 = t0 + dt;
t2 = t1 + dt;

init.t0 = t0;
init.t1 = t1;
init.t2 = t2;
init.u0 = structts2vec(s.mpcsoln, cv.unames, t0);
init.u1 = structts2vec(s.mpcsoln, cv.unames, t1);
init.x1 = structts2vec(s.mpcsoln, cv.xnames, t1);
init.x2 = structts2vec(s.mpcsoln, cv.xnames, t2);
init.e1 = structts2vec(s.mpcsoln, cv.enames, t1);
init.e2 = structts2vec(s.mpcsoln, cv.enames, t2);
init.pcurrt0 = structts2vec(s.mpcsoln, {'pcurrt'}, t0);

eqfds = {'psizr','psiapp','psipla','iv','ic','pcurrt','psibry'};
init.eq1 = structts2struct(s.mpcsoln, eqfds, t1);
init.eq2 = structts2struct(s.mpcsoln, eqfds, t2);

init.eq1 = reshape_fds(init.eq1, tok);
init.eq2 = reshape_fds(init.eq2, tok);

end

function eq = reshape_fds(eq, tok)

  eq.pcurrt = reshape(eq.pcurrt, tok.nz, tok.nr);
  eq.psizr  = reshape(eq.psizr,  tok.nz, tok.nr);
  eq.psiapp = reshape(eq.psiapp, tok.nz, tok.nr);
  eq.psipla = reshape(eq.psipla, tok.nz, tok.nr);
end