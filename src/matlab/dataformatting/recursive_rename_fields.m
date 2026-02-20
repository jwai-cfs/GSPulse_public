function obj = recursive_rename_fields(obj, fname_old, fname_new)
  obj = recursive_struct_cell_mods(obj, @rename_fields, fname_old, fname_new);
end

function s = rename_fields(s, fname_old, fname_new)
  if isstruct(s)
    fnames = fieldnames(s);
    for i = 1:length(fnames)
      fname = fnames{i};
      if strcmp(fname, fname_old)
        s.(fname_new) = s.(fname_old);
        s = rmfield(s, fname_old);
      end
    end
  end
end