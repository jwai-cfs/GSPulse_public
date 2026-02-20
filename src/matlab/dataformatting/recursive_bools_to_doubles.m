function obj = recursive_bools_to_doubles(obj)
  obj = recursive_struct_cell_mods(obj, @bool_to_double);
end

function val = bool_to_double(val)
  if islogical(val)
    val = double(val);
  end
end