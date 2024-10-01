function vec = structts2vec(s, fdnames, time)
% Helper function for formatting a struct of timeseries (or
% timeseries-like) into a vector. The inverse of this procedure is
% vec2structts.m
% 
% Example: 
% t = linspace(0,2*pi, 500)';
% s.sig1 = timeseries(sin(t), t);
% s.sig2 = timeseries(cos(t), t);
% fdnames = {'sig1', 'sig2'};
% teval = [0 pi/2]; 
% vec = structts2vec(s, fdnames, teval);  % evaluate it at t=0,pi/2 
% vec  % (should be [0 1 1 0])

x = cell(length(fdnames),1);
for i = 1:length(fdnames)  
  fd = fdnames{i};

  % this is an edge case where Time and Data only contain a single value
  % and the normal interp1 behavior fails
  if (length(s.(fd).Time)) == 1 && (s.(fd).Time == time)
    x{i} = s.(fd).Data;
  else
    x{i} = interp1(s.(fd).Time, s.(fd).Data, time(:), 'linear', 'extrap'); 
  end


end

x = horzcat(x{:})';

vec = x(:);









