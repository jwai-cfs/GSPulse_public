function plot_lim(tok, varargin)
% =========================================================================
% Description: 
%   plot the tokamak limiter
%
% Inputs: 
%  tok - tokamak geometry object
%  varargin - plotting specifications
%
% Outputs:
%   None
% =========================================================================

if isempty(varargin)
  varargin = {'color', [0.4 0.4 0.4], 'linewidth', 1.5};
end

[rl, zl] = close_curve(tok.rl, tok.zl);

hold on
plot(rl, zl, varargin{:});

axis equal
axis([min(tok.rg) max(tok.rg) min(tok.zg) max(tok.zg)])








