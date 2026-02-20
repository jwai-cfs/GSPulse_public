function A = copyfields(A, B, fields2copy, overwrite)
% Description: 
% Copies the fields of B to the fields of A. 
%
% fields2copy should be a cell of strings corresponding to names of fields
% in B to copy to A. Use fields2copy={} to copy all fields of B
%
% If overwrite=1, the fields of B will overwrite contents of corresponding 
% field in A. If overwrite=0, only fields of B that do not exist in A will
% be copied. 

% Inputs: 
%  A - base struct
%  B - struct from which to copy data
%  fds2copy -fields to copy from B, use fds2copy={} to copy all fields of B
%  overwrite - if overwrite=1, the fields of B will overwrite contents of 
%   corresponding field in A. If overwrite=0, only fields of B that do not
%   exist in A will be copied.
%
% Output: 
%  A - struct with the requested fields from B copied over

if nargin < 4, overwrite = false; end

if isempty(fields2copy), fields2copy = fieldnames(B); end

for i = 1:length(fields2copy)
  fn = fields2copy{i};  
  if ~isfield(A,fn) || overwrite
    A.(fn) = B.(fn);
  end
end
