function s = rmfields(s, fds)
% Description: 
% remove multiple fields at a time from a struct 
%
% Inputs: s - input struct, fds - fields to be removed
% Outputs: s - struct with fds purged

for i = 1:length(fds)
  fd = fds{i};
  if isfield(s, fd)
    s = rmfield(s, fd);
  end
end