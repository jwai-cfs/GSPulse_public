function [fighandle, axs] = plot_structts(s, fdnames, fighandle, nrows, ncols, varargin)
% Description: 
% plot data stored that is stored as a struct of timeseries (or
% timeseries-like, with Time and Data fields). 
%
% Example: 
% t = linspace(0,2*pi)';
% s.sig1 = timeseries(sin(t), t);
% s.sig2 = timeseries(cos(t), t);
% fdnames = {'sig1', 'sig2'};
% h = plot_structts(s, fdnames, 1, [], 'linewidth', 2); 

N = length(fdnames);
if ~exist('nrows','var') || isempty(nrows), nrows = min(3,N); end
if ~exist('ncols','var') || isempty(ncols), ncols = ceil(N/nrows); end
if ~exist('fighandle','var') || isempty(fighandle), fighandle = figure; end


figure(fighandle);
axs = [];

for i = 1:N
  fd = fdnames{i};
  axs(i) = subplot(nrows, ncols, i);
  hold on
  grid on
  plot( s.(fd).Time, s.(fd).Data, varargin{:});
  title(fd, 'fontsize', 16, 'fontweight', 'bold')
  if isfield(s.(fd), 'Units')
    ylabel(s.(fd).Units, 'fontsize', 13)
  end
  xlabel('Time [s]', 'fontsize', 13)
end
linkaxes(axs, 'x')
