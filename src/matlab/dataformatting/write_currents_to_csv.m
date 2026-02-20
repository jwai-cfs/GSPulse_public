function write_currents_to_csv(fn, signals, ccnames)

t = signals.ic.Time;
ic = signals.ic.Data;
colheaders = [{'Time [sec]'}; strcat(ccnames(:), 'current [A/turn]')];
T = array2table([t ic], 'VariableNames', colheaders);

f = fileparts(fn);
if ~isfolder(f)
  mkdir(f); 
end
writetable(T, fn);
