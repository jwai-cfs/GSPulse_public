function s = smooth_structts(s, dt_window)

% smooth each timeseries
fds = fieldnames(s);
for i = 1:length(fds)
   fd = fds{i};
   s.(fd) = smooth_ts(s.(fd), dt_window);
end