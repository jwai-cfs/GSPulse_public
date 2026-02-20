function s = variables2struct(varargin)
% Description: creates a struct with each field and fieldname given by the input args
%
% Example:
% a = 1;
% b = 2;
% s = variables2struct(a,b);

s = struct;
for i = 1:nargin
  s.(inputname(i)) = varargin{i};
end
end