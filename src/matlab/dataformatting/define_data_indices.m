function cv = define_data_indices(optimization_signals, tok)
% Description: 
% define data indices for the (c)ontrolled (v)ariables
%
% Inputs: 
%
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
cv.y_names = {};
cv.u_name = {};
cv.e_names = {};
idx = 0;


for i = 1:length(optimization_signals)
  
  sig = optimization_signals{i};

  % define data indices for outputs, y, and errors, e
  if ~strcmp(sig.calc_type, 'voltage')
    n = size(sig.wt.Data, 2);
    if n > 0
      idx = idx(end) + (1:n);
  
      cv.iy.(sig.name) = idx;
      cv.y_names{end+1} = sig.name;
      % cv.y_descriptions.(sig.name) = sig.description; 
  
      cv.ie.([sig.name '_err']) = idx;
      cv.e_names{end+1} = [sig.name '_err'];
      cv.e_descriptions.(sig.name) = ['Error (target - value) for signal ' sig.name];
    end
  % data indices for inputs, u
  else
    cv.u_name{end+1} = sig.name;
    cv.iu.(sig.name) = 1:size(sig.wt.Data,2);
    % cv.u_description.(sig.name) = sig.description;
  end
end


% data indices for state variables, x
cv.x_names = {'ic', 'ivb'};
cv.ix.ic = 1:tok.nc;
cv.ix.ivb = tok.nc + (1:tok.nv);
cv.x_descriptions.ic = 'Coil currents.';
cv.x_descriptions.ivb = 'Vessel currents (compressed format)';
