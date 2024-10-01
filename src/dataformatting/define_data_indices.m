 function cv = define_data_indices(settings, targs, tok)
% Description: 
% define data indices for the (c)ontrolled (v)ariables
%
% Inputs: 
%  settings - settings object, see help _define_settings.m
%  targs    - targs object, see help _define_targets.m
%  tok      - tok object, see help _define_tok.m
%
% Outputs:
%  cv - struct with information on the data indices of the controlled
%  variables
% 
% Additional info:
%  The controlled variables can generally be grouped into outputs (y), 
%  states (x), and inputs (u). There can be multiple subcategories for
%  each. cv.iy, cv.ix, and cv.iu define databuses for organizing and
%  referencing these controlled variables by name and index. 

cv = struct;


% the outputs y are defined by fds2control
cv.ynames = settings.fds2control;
idx = 0;
for k = 1:length(cv.ynames)    
  varname = cv.ynames{k};
  n = size(targs.(varname).Data, 2);  
  idx = idx(end)+1:idx(end)+n;  
  cv.iy.(varname) = idx;
end 
for i = 1:tok.nc
  cv.iy.(tok.ccnames{i}) = cv.iy.ic(i);  
end


% the errs, e (same data indices as y)
cv.enames = strcat(settings.fds2control, '_err');
iy_fds = fieldnames(cv.iy);
ie_fds = strcat(iy_fds, '_err');
for i = 1:length(ie_fds)
  cv.ie.(ie_fds{i}) = cv.iy.(iy_fds{i});
end


% the states x are the coil currents and vessel current modes
cv.xnames = {'ic', 'ivb'};
cv.ix.ic = 1:tok.nc;
cv.ix.ivb = tok.nc + (1:settings.nvessmodes);
for i = 1:tok.nc
  cv.ix.(tok.ccnames{i}) = i;  
end
cv.xdesc.ic = 'Coil currents.';
cv.xdesc.ivb = 'Vessel currents (compressed format)';


% the inputs are the voltages on the active coils (subset of the coils)
cv.unames = {'v'};
for i = 1:tok.nc
  nam = [tok.ccnames{i} '_V'];
  cv.iu.(nam) = i;
end
cv.iu.v = 1:tok.nc;
cv.udesc.v = 'Voltage in the active coil circuits.';