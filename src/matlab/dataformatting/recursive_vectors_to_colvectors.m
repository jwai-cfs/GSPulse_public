function obj = recursive_vectors_to_colvectors(obj)
  obj = recursive_struct_cell_mods(obj, @vectors_to_colvectors);
end

function val = vectors_to_colvectors(val)
  if isvector(val) && ~ischar(val) && ~isa(val, 'function_handle')
    val = val(:);
  end
end