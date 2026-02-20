function vec = struct2vec(s, fdnames)
% Description: 
% convert a struct of data fields into a single vector (also works with 2D
% data to output a matrix). 
%
% Example:
% s.a = [1 1 1 1];
% s.b = [2 2 2 2];
% s.c = [3 3 3 3; 4 4 4 4];
% fdnames = {'a','b','c'};
% vec = struct2vec(s, fdnames)

x = cell(length(fdnames),1);
for i = 1:length(fdnames)    
  x{i} = s.(fdnames{i});
end
vec = vertcat(x{:});



