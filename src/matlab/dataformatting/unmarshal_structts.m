function s = unmarshal_structts(s, bus_name, signal_names)

s = check_structts_dims(s);

assert(length(signal_names) == size(s.(bus_name).Data,2), 'Data size inconsistent')

for i = 1:length(signal_names)
  s.(signal_names{i}).Time = s.(bus_name).Time;
  s.(signal_names{i}).Data = s.(bus_name).Data(:,i);
end
