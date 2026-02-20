function ts = blend_ts(ts1, ts2, alpha)
% Description:
%  blends two timseries together according to the waveform specified in
%  alpha
%
% Inputs:
%  ts1 - input timeseries 1
%  ts2 - input timeseries 2
%  alpha - time dependent blend fraction
%
% Outputs:
%  blended timeseries
%  
% Example:
%  ts1.Time = 1:10;
%  ts1.Data = cos(1:10)';
%  ts2.Time = 11:20;
%  ts2.Data = cos(11:20)';
%  ts = blend_ts(ts1, ts2, alpha);
t = union(ts1.Time(:), ts2.Time(:));
t = sort(t);
a = interp1hold(alpha.Time, alpha.Data, t);

data1 = interp1hold(ts1.Time, ts1.Data, t);
data2 = interp1hold(ts2.Time, ts2.Data, t);
data1(isnan(data1)) = 0;
data2(isnan(data2)) = 0;
data1 = sparse(diag(a)) * data1;
data2 = sparse(diag(1-a)) * data2;

ts.Time = t;
ts.Data = data1 + data2;



