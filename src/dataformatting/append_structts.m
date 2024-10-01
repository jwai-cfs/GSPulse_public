function sglobal = append_structts(varargin)
% Description:
%  Concatenates the timeseries data from multiple struccts (struct of
%  timeseries) objects. 
%
% Inputs:
%  variable of number of struccts objects  
%
% Outputs:
%  concatenated struct of timeseries
%  
% Example:
%  s1.a.Time = 1:10;
%  s1.a.Data = cos(1:10)';
%  s2.a.Time = 11:20;
%  s2.a.Data = cos(11:20)';
%  s = append_structts(s1, s2);
%
% Restrictions:
%  does not do any sorting, does not check for consistency of variable
%  names in the structts objects, does not verify data dimensions

sglobal = varargin{1};

for i = 2:nargin
  s = varargin{i}; 
  fds = fieldnames(s);
  for j = 1:length(fds)
    fd = fds{j};    
    sglobal.(fd).Time = [sglobal.(fd).Time(:); s.(fd).Time(:)];
    sglobal.(fd).Data = [sglobal.(fd).Data; s.(fd).Data];        
  end
end
