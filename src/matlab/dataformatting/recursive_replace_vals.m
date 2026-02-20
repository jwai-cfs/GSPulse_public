function obj = recursive_replace_vals(obj, check_val, replace_val)
  obj = recursive_struct_cell_mods(obj, @replace_val, check_val, replace_val);
end

function val = replace_val(val, check_val, replace_val)
  if isequaln(val, check_val)
    val = replace_val;
  end  
end