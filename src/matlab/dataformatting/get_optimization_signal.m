function [sig, i] = get_optimization_signal(optimization_signals, name)

sig = struct;
n = length(optimization_signals);
for i = 1:n
  if strcmp(optimization_signals{i}.name, name)
    sig = optimization_signals{i};
    return
  end
end

i = nan;
warning("Could not find signal '%s' in list of optimization signals", name)
