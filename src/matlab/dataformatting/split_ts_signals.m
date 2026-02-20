function ts1 = split_ts_signals(ts, sigs)
% Description: 
%  split a timeseries containing matrix data into a timeseries containing
%  multiple signals of 1-D vector data
% 
% Inputs: ts - timeseries with matrix data, sigs - names of the signals to
% split the matrix data into
%
% Outputs: ts - timeseries with multiple signals each containing vector data

ts1 = struct;
sigs = cellstr(sigs);
ts = check_ts_dims(ts);

for i = 1:length(sigs)
  ts1.(sigs{i}).Time = ts.Time;
  ts1.(sigs{i}).Data = ts.Data(:,i);
end
