function [wts, dwts, d2wts, targs, islocked] = assemble_weights_and_targets(optimization_signals)

wts = struct;
dwts = struct;
d2wts = struct;
targs = struct;
islocked = struct; 

for i = 1:length(optimization_signals)
  
  sig = optimization_signals{i};

  wts.(sig.name) = sig.wt;
  dwts.(sig.name) = sig.dwt;
  d2wts.(sig.name) = sig.d2wt;
  targs.(sig.name) = sig.targ;    
  
  if isfield(sig, 'constrain')
    islocked.(sig.name) = sig.constrain;
  else
    islocked.(sig.name).Time = sig.targ.Time;
    islocked.(sig.name).Data = zeros(size(sig.targ.Data));
  end

end
