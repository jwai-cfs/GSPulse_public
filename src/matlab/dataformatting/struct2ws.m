function struct2ws(s)
% Description: 
%   makes all the fields of s available as workspace variables (or function
%   variables, if called from within a function)
% 
% Inputs: s - struct with data fields
% Outputs: fields of s are now available as variables

nams = fieldnames(s);
for i = 1:length(nams)
  nam = nams{i};
  assignin('caller', nam, s.(nam));
end
end