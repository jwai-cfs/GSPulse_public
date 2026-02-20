function obj = recursive_struct_cell_mods(obj, modfun, varargin)
% =========================================================================
% Recursively operate on the {fields of a struct, or entries within a cell
% array} by applying the @modfun to each sub-entry. 
%
% Inputs: obj - the object (cell array or struct) to be modified in place
%         modfun - function handle that operates on the struct field values
%                  or cell values
%         varargin - additional arguments passed to modfun
%
% Outputs: obj - modified obj
% =========================================================================
obj = modfun(obj, varargin{:});  % apply modifications to root obj 

if iscell(obj)
  
  % loop through the entries of a cell array
  for i = 1:length(obj)
    sub_obj = obj{i}; 
  
    % recursion step 
    if isstruct(sub_obj) || iscell(sub_obj)
      obj{i} = recursive_struct_cell_mods(sub_obj, modfun, varargin{:});
    
    % apply modifications to sub-object
    else
      obj{i} = modfun(sub_obj, varargin{:});
    end
  end

elseif isstruct(obj)

  % loop through the fields of a struct
  fnames = fieldnames(obj);
  for i = 1:length(fnames)
    fname = fnames{i};
    sub_obj = obj.(fname);
    
    % recursion step
    if isstruct(sub_obj) || iscell(sub_obj)
      obj.(fname) = recursive_struct_cell_mods(sub_obj, modfun, varargin{:});
    
    % apply modifications to sub-object
    else      
      obj.(fname) = modfun(sub_obj, varargin{:});
    end
  end
end