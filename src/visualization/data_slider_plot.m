function fig = data_slider_plot(ts, step, x)
% =========================================================================
% Description: 
%  create a slider plot for a timeseries that contains 1-D data vs time
% (e.g. profile data vs time)
% 
% Inputs: 
%  ts - timeseries with 1-D data vs time
%  step - int, step size for the slider bar
%  x - optional, x-axis for the plot
%
% Outputs: 
%  fig - figure handle
% =========================================================================

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
s.Callback = {@sliderCallback, ts, x};

plot_(s.Value, ts, x)

% slider callback
function sliderCallback(src, event, ts, x)
  plot_(round(src.Value), ts, x)
end


% plot shape targets
function plot_(i, ts, x)
    
  t = ts.Time(i);
  y = squeeze(interp1(ts.Time, ts.Data, t));
  if isempty(x), x = 1:length(y); end

  cla
  hold on
  plot(x,y)
  title(t)
  grid on
  ylim([min(ts.Data(:)) max(ts.Data(:))])
  drawnow  

end
end