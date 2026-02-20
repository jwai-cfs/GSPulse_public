function s = interp1_structts(s, tnew, varargin)

% interpolate each timeseries
fds = fieldnames(s);
for i = 1:length(fds)
   ts = s.(fds{i});
   ts2.Time = tnew;
   ts2.Data = interp1(ts.Time, ts.Data, tnew, varargin{:});
   s.(fds{i}) = ts2;
end

