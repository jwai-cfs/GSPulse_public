function s = validate_structts_timing(s, required_tspan)

s = check_structts_dims(s);
signal_names = fieldnames(s);
assert(numel(required_tspan) == 2)

for i = 1:length(signal_names)
  signal_name = signal_names{i};  

  if isfield(s.(signal_name), 'Time') && isfield(s.(signal_name), 'Data')
    % Check that time is sorted
    assert(issorted(s.(signal_name).Time), ['Time for signal "%s" is not' ...
      ' sorted'], signal_name);
    
    % Check that time points are unique
    assert(numel(s.(signal_name).Time) == numel(unique(s.(signal_name).Time)), ...
      'Time points in signal "%s" are not unique', signal_name);

    actual_tspan = s.(signal_name).Time([1 end]);
    actual_tarray = s.(signal_name).Time;
  
    % Check that signal time spans tspan
    assert(actual_tspan(1) <= required_tspan(1) && actual_tspan(2) >= required_tspan(2), ...
      ['Signal "%s" is specified with time range [%.3f, %.3f] which does ' ...
      'not span the required time interval [%.3f, %.3f].'], ...
      signal_name, actual_tspan(1), actual_tspan(2), required_tspan(1), required_tspan(2));
  

    % Check dimension match of Time and Data
    assert(size(s.(signal_name).Time,1) == size(s.(signal_name).Data,1), ...
      'Time and Data dimensions inconsistent for signal "%s" (Time=%d, Data=%d)', ...
       signal_name, size(s.(signal_name).Time,1), size(s.(signal_name).Data,1));
  
    
    % Handle infs (first and last time points allowed to be inf)
    % ...........................................................

    isneginf = @(x) isinf(x) && x < 0; 
    isposinf = @(x) isinf(x) && x > 0;

    if isneginf(actual_tspan(1)) 

      % data window that uses the inf is not used, remove that data
      if actual_tarray(2) <= required_tspan(1)
        s.(signal_name).Time(1) = [];
        s.(signal_name).Data(1,:) = [];
      
      % data window with inf is used, perform checks and re-map the inf to
      % the required_tspan
      else
        assert( all(isequaln(s.(signal_name).Data(1,:), s.(signal_name).Data(2,:))), ...
        ['Signal "%s" time range starts with -inf, but interpolation ' ...
        'cannot be performed because adjacent data value is not equal.'], signal_name);
    
        s.(signal_name).Time(1) = required_tspan(1);
      end
    end

    if isposinf(actual_tspan(end)) 

      % data window that uses the inf is not used, remove that data
      if actual_tarray(end-1) >= required_tspan(end) 
        s.(signal_name).Time(end) = [];
        s.(signal_name).Data(end,:) = [];
      
      % data window with inf is used, perform checks and re-map the inf to
      % the required_tspan
      else
        assert( all(isequaln(s.(signal_name).Data(end,:), s.(signal_name).Data(end-1,:))), ...
        ['Signal "%s" time range ends with inf, but interpolation ' ...
        'cannot be performed because adjacent data value is not equal.'], signal_name);

        s.(signal_name).Time(end) = required_tspan(end);
      end
    end
  end
end
