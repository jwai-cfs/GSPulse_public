

function fig = data_slider_plot(ts, ts2, step, x)

if ~exist('x', 'var'), x = []; end
if ~exist('step', 'var'), step = 1; end

% ensure uniform spacing
t_uniform = linspace(ts.Time(1), ts.Time(end), length(ts.Time));
ts.Data = interp1(ts.Time, ts.Data, t_uniform);
ts.Time = t_uniform;

% create figure
fig = figure;
fig.Position = [749 178 662 499];
ax = axes(fig);
ax.Box = 'on';
ax.Position = [0.13 0.2 0.775 0.7];

% slider button
s = uicontrol(fig, 'style', 'slider');
s.Units = 'normalized';
s.Position = [0.15 0.05 0.6 0.05];
s.Min = 1;
s.Max = length(ts.Time);
s.Value = 1;
s.SliderStep = step * [1 1] / (s.Max - s.Min);
s.Callback = {@sliderCallback, ts, ts2, x};

plot_(s.Value, ts, ts2, x)

% slider callback
function sliderCallback(src, event, ts, ts2, x)
  plot_(round(src.Value), ts, ts2, x)
end


% plot shape targets
function plot_(i, ts, ts2, x)
    
  t = ts.Time(i);
  y = squeeze(interp1(ts.Time, ts.Data, t));
  y2 = squeeze(interp1(ts2.Time, ts2.Data, t));
  if isempty(x), x = 1:length(y); end

  cla
  hold on
  plot(x,y)
  plot(x,y2,'--')
  title(t)
  grid on
  ylim([min(ts.Data(:)) max(ts.Data(:))])
  drawnow  

end
end