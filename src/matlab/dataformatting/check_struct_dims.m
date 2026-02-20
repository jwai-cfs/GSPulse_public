function s = check_struct_dims(s, mode)
% Description:
% for each field of the struct s, if the field corresponds to a numeric
% vector, make the vector be a column vector (default) or row vector
%
% Inputs:
%  s - input data struct
%  mode - 'cols' or 'rows', which way to force 1-D vectors
%
% Outputs:
%  s - struct, any 1-D numeric vectors have been transposed to be either
%  columns or rows according to 'mode'
if ~exist('mode', 'var'), mode = 'cols'; end

fds = fieldnames(s);
for i = 1:length(fds)
  fd = fds{i};
  if isvector(s.(fd))    
    switch mode
      case 'cols'
        s.(fd) = s.(fd)(:);        
      case 'rows'
        s.(fd) = s.(fd)(:)';
      otherwise
        error('Unknown mode')
    end
  end  
end