function obj = recursive_vals_to_doubles(obj)
  obj = recursive_struct_cell_mods(obj, @val_to_double);
end

function val = val_to_double(val)
  if isnumeric(val)
    val = double(val);
  end
end