function obj = recursive_remove_nan_and_empty_fields(obj)
  obj = recursive_struct_cell_mods(obj, @remove_nan_and_empty_fields);
end

function val = remove_nan_and_empty_fields(val)   
  if isstruct(val)
    fds = fieldnames(val);
    for i = 1:length(fds)
      subval = val.(fds{i});
      if isnumeric(subval) && numel(subval) <= 1 
        if isempty(subval) || isnan(subval)
          val = rmfield(val, fds{i});
        end
      end        
    end
  end 
end