function time_data = parse_time_data(interval_t)

time_data = struct;
time_data.n_intervals = length(interval_t);

for i = 1:time_data.n_intervals
  t1N = interval_t{i};       % time array from 1 to N
  t1N = t1N(:);              
  N = length(t1N);           % number of steps
  t1 = t1N(1);               % time at timestep 1
  tN = t1N(end);             % time at timestep N  

  if N > 1
    dt1N = diff(t1N);           % difference in time steps for steps 1 to N
    t0 = t1 - dt1N(1);          % time at timestep 0, inferred
    tNm1 = t1N(end-1);          % time at timestep N-1
    t0Nm1 = [t0; t1N(1:end-1)]; % time array from 0 to N-1
  else
    dt1N = [];
    t0 = t1;
    tNm1 = [];
    t0Nm1 = t0;
  end

  % Check that time is sorted and unique
  assert(issorted(t1N), 'Global solution time interval_t{%d} is not sorted', i);
  assert(numel(t1N) == numel(unique(t1N)), ['Global solution time ' ...
    'interval_t{%d} time points are not unique'], i); 

  time_data.interval(i) = variables2struct(t1N, N, dt1N, t0, t0Nm1, t1, tN, tNm1);


  % Validate timing
  if i > 1
    if length(t1N) >= 2
      if t1N(2) > time_data.interval(i-1).tN
        str = ['\nTime intervals in settings.interval_t do not overlap by 2 or more samples. \n' ...
          'Interval %d ends at %.3f sec and interval %d begins with ' ...
          '[%.3f, %.3f] sec. \nIncrease overlap time window. \n'];

        error(str, i-1, time_data.interval(i-1).tN, i, t1N(1), t1N(2))
      end
    end
  end
end
