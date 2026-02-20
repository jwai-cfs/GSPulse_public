function s = blend_structts(s1, s2, alpha, fds)
% Description:
%  blends the waveforms of data stored in the s1 and s2 structts objects,
%  according to the blend fraction specified by the alpha waveform
%
% Inputs:
%  s1 - struct of timeseries object
%  s2 - struct of timeseries object
%  alpha - blend fraction timeseries object
%  fds - [optional], which fields to blend
%
% Outputs:
%  blend struct of timeseries
%  
% Example:
%  s1.a.Time = 1:100;
%  s1.a.Data = cos(1:100)';
%  s2.a.Time = 1:100;
%  s2.a.Data = sin(1:100)';
%  alpha.Time = 1:100;
%  alpha.Data = linspace(1,0,100)';
%  s = blend_structts(s1, s2, alpha, {'a'});
%
% Restrictions:
%  does not do any sorting, does not check for consistency of variable
%  names in the structts objects, does not verify data dimensions

if ~exist('fds', 'var') || isempty(fds)
  fds = intersect(fieldnames(s1), fieldnames(s2));
end

if isnumeric(alpha)
  a = alpha;
  alpha = struct;
  alpha.Time = [0 1];
  alpha.Data = [a a];
end

s = struct; 
for i = 1:length(fds)
  fd = fds{i};
  s.(fd) = blend_ts(s1.(fd), s2.(fd), alpha);
end