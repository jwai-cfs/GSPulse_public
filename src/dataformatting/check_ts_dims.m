function ts = check_ts_dims(ts)
% Description:
%
% transpose timeseries data such that the first dimension of the data has
% the same length as the time
%
% ts is a struct with fields that have Data and Time subfields
%    ts.Data = ...
%    ts.Time = ...
%
% This function verifies that ts.Time is a column vector and that the
% first dimension of ts.Data corresponds to time, and transposes any
% results that dont fit this pattern. 
% 
% Inputs:
%   ts - struct of timeseries data
% 
% Outputs:
%   ts - struct with the timeseries data formatted in a consistent manner
%
% Restrictions:
% does not work for data with 3 or more dimensions. 
if isfield(ts, 'Time') && isfield(ts, 'Data')
  ts.Time = ts.Time(:);
  n = length(ts.Time);
  sz = size(ts.Data);

  if sz(1) ~= n
    if sz(2) == n && length(sz) == 2
      ts.Data = ts.Data';
    end
  end
end