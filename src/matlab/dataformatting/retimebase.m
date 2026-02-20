function s = retimebase(s, t)
% Description: 
% interpolate a struct of timeseries (or timeseries-like) onto new
% timebase, t
%
% Inputs: 
%    s - input struct of timeseries
%    t - new timebase to interpolate onto
%
% Outputs: 
%    s - struct of timeseries with new timebase
%
% Example:
% t = linspace(0,2*pi,100)';
% s.sig1 = timeseries(sin(t), t);
% s.sig2 = timeseries(cos(t), t);
% tnew = t(1:5:end);
% s = retimebase(s, t)
%
% Warnings: this uses interp1hold, which can quietly extrapolate data

% interpolate onto new timebase
fds = fieldnames(s);
for i = 1:length(fds)
  fd = fds{i};
  if isstruct(s.(fd)) && isfield(s.(fd), 'Data')
    s.(fd).Data = interp1hold(s.(fd).Time, s.(fd).Data, t);
    s.(fd).Time = t;
  elseif isa(s.(fd), 'timeseries')
    s.(fd) = resample(s.(fd), t);
  end

end
