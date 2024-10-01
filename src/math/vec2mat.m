function A = vec2mat(vec, nrows, shift)
% =========================================================================
% Description:
% Creates a matrix by repeating a row vector on subsequent rows of the
% matrix. 
%
% Inputs: 
%   vec - row vector to repeat
%   nrows - number of rows in the output matrix
%   shift - the number of elements to shift the vector by in each
%           subsequent row
% 
% Example: this example creates a finite-difference Laplacian matrix
% vec = [1 -2 1];
% shift = 1;
% nrows = 4;
% A = vec2mat(vec, nrows, shift)
% =========================================================================
vec = vec(:)';
l = length(vec);
A = zeros(nrows, nrows + l*(shift-1));

for i = 1:nrows
  j = (i-1)*shift + (1:l);
  A(i,j) = vec;
end
