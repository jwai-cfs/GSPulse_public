function obj = recursive_struct_array_to_cell(obj)
  obj = recursive_struct_cell_mods(obj, @struct_array_to_cell);
end

function val = struct_array_to_cell(val)
  
  if isstruct(val) && numel(val) > 1
    newval = cell(size(val));
    for i = 1:numel(val)
      newval{i} = val(i);
    end
    val = newval;
  end
  
end