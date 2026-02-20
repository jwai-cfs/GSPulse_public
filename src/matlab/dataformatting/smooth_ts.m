function ts2 = smooth_ts(ts, dt_window)

dt = mean(diff(ts.Time));
nsm = ceil(dt_window/dt);

ts2.Time = ts.Time;
ts2.Data = smoothdata(ts.Data, 1, 'gaussian', nsm);