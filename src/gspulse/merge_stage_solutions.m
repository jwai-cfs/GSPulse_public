function [tglobal, eqs, signals] = merge_stage_solutions(stagesoln, timedata)
% =========================================================================
% Description: 
%  merge the stage solutions, with overlapping time intervals, into one
%  consistent global solution with a single time interval
%
% Inputs: 
%  stagesoln - solution struct from GSPulse for that particular stage
%  timedata - timing information, see extract_timedata.m
% 
% Outputs: 
%  tglobal - vector, global timebase
%  eqs - cell array, containing equilibria for each time in tglobal
%  signals - all of the MPC optimization waveforms, mapped onto tglobal.
%
% =========================================================================
calc_signals = isfield(stagesoln{1}, 'mpcsoln');

nint = timedata.n_intervals;

tglobal = [];
eqs = {};

for kint = 1:nint     
  
  if kint < nint
    tnext = timedata.interval(kint+1).t;
  else
    tnext = inf(3,1);
  end
  ikeep = timedata.interval(kint).t + sqrt(eps) < tnext(3);
  if kint > 1
    ikeep(1:2) = 0;
  end

  t = timedata.interval(kint).t(ikeep);
  tglobal = [tglobal; t];
  eqs = [eqs; stagesoln{kint}.eqs(ikeep)];
  
  if calc_signals
    mpcsolns{kint} = retimebase(stagesoln{kint}.mpcsoln, t);
    mpcsolns{kint} = check_structts_dims(mpcsolns{kint});
  end
end

if calc_signals
  signals = append_structts(mpcsolns{:});
else
  signals = [];
end