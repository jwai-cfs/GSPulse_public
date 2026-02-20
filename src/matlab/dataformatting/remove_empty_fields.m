function s = remove_empty_fields(s)

fds = fieldnames(s);
for i = 1:length(fds)
  if isempty(s.(fds{i}))
    s = rmfield(s, fds{i});
  end
end
