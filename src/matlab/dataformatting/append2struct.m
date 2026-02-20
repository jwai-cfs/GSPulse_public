function s = append2struct(s, varargin)  
% DESCRIPTION: append variables to a struct as fields of that struct
%
% EXAMPLE:
% s.a = 1; 
% b = 2; 
% s = append2struct(s, b);  % returns struct with s.a=1, s.b=2

  for i = 2:nargin
    s.(inputname(i)) = varargin{i-1};
  end
end