function s = check_structts_dims(s)
% Description:
%
% Ensure consistent formatting of a struct of timeseries object. Transposes
% each data field in each timeseries such that the first dimension of 
% each data field has the same length as the time field. 
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
%
% Example: 
% t = linspace(0,1);
% s.a.Time = t;        % row vector
% s.a.Data = sin(t)';  % col vector
% s.b.Time = t';       % col vector 
% s.b.Data = cos(t);   % row vector
%
% s = check_structts_dims(s);
% s.a    % all Time and Data fields are now col vectors
% s.b  

fds = fieldnames(s);
for i = 1:length(fds)
  fd = fds{i};
  s.(fd) = check_ts_dims(s.(fd));
end

      
        
        
      
        






