function obj = recursive_char_to_cellstr(obj)
  obj = recursive_struct_cell_mods(obj, @char_to_cellstr);
end

function val = char_to_cellstr(val)
  if ischar(val)
    if isvector(val)
      if iscolumn(val)
        val = val';       
      end
    else
      val = cellstr(val);  
    end
  end
end

