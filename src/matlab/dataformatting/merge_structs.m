function s = merge_structs(varargin)
% Description:
% Copies all the fields of each input structure to a new structure. 
%
% Example:
% s1.a = 1
% s2.b = 2
% s = merge_structs(s1, s2)
% 
% Warnings: note that this can overwrite data

s = struct;
for i = 1:nargin
  tmp = varargin{i}; 
  fds = fieldnames(tmp);
  for j = 1:length(fds)
    fd = fds{j};
    s.(fd) = tmp.(fd);
  end
end









