function varargout = nan_or_empty_to_zero(varargin)
% Description: 
%  returns the same as the inputs, except if the inputs contain nan or empty
%  values, then returns 0
%
% Example:
% a = [];
% b = [2 0];
% c = [nan nan nan];
% [a, b, c] = nan_or_empty_to_zero(a, b, c);
% returns a = 0, b = [2 0]; c = [0 0 0];

  varargout = cell(1,nargin);

  for i = 1:nargin
    x = varargin{i};
    if all(isnan(x))
      x = zeros(size(x));
    end
    if isempty(x)
      x = 0;
    end
    varargout{i} = x;
  end  

end


