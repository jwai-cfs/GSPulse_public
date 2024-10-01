function timedata = extract_timedata(interval_t)
% =========================================================================
% Description: 
%  Defines a bunch of auxiliary timing data based on the time intervals
%   specified by settings.interval_t
%
% Inputs: 
%  interval_t - cellarray with multiple timeintervals supplied, see example
%
% Outputs: 
%  timedata - struct array with more auxiliary timing data defined (such as
%  the dt spacing, first time, end time, etc)
%
% Restrictions: 
% each interval must have uniform spacing
%
% Example: 
% dt = 0.2;
% interval_t{1} = 0.2:dt:9;
% interval_t{2} = 5:dt:9;
% interval_t{3} = 9:dt:12;
% timedata = extract_timedata(interval_t)
% =========================================================================
timedata = struct;
n_intervals = length(interval_t);
timedata.n_intervals = n_intervals;

% timing data for each individual interval
for kint = 1:n_intervals

  t1N = interval_t{kint};                 % timebase from 0 to N
  t1N = t1N(:);
  if length(t1N) > 1
    dt = mean(diff(t1N));                 % discretization time
  else
    dt = 0;
  end
  t1 = t1N(1);                            % time at index 1
  t0 = t1 - dt;                           % time at index 0
  tN  = t1N(end);                         % time at index N;
  t = t1N;                                % the most commonly used timebase, shorthand for t1N
  t0Nm1 = t1N - dt;                       % timebase from 0 to N-1
  N = length(t);

  timedata.interval(kint) = ...
    variables2struct(t1N, dt, t0, t1, tN, t1N, t0Nm1, t, N);

  if kint > 1
    if length(t) >= 2
      if t(2) > timedata.interval(kint-1).tN
        str = ['\nTime intervals in settings.interval_t do not overlap by 2 or more samples. \n' ...
               'Interval %d ends at %.3f sec and interval %d begins with ' ...
               '[%.3f, %.3f] sec. \nRecommend to increase overlap (see define_general_settings.m) \n'];

        warning(str, kint-1, timedata.interval(kint-1).tN, kint, t(1), t(2))
      end
    end
  end
end