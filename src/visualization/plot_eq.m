function plot_eq(eq, tok, varargin)
% =========================================================================
% Description: 
%   makes a simple equilibrium plot of the boundary contour
%
% Inputs: 
%  eq - GSpulse equilibrium struct
%  tok - tokamak geometry struct
%  varargin - plotting specifications for the call to contour
%
% Outputs: 
%  None
% =========================================================================

if isfield(eq, 'psizr')
  psizr = eq.psizr;
elseif isfield(eq, 'Fx')
  psizr = eq.Fx;
else
  warning('Could not read psizr')
end

if isfield(eq, 'psibry')
  psibry = eq.psibry;
elseif isfield(eq, 'FB')
  psibry = eq.FB;
else
  warning('Could not read psibry')
end

if isfield(tok, 'rg')
  rg = tok.rg;
  zg = tok.zg;
else
  warning('could not read (rg,zg)')
end

hold on
plot(tok.rl, tok.zl, 'k')
contour(rg,zg,psizr,[psibry psibry], varargin{:});
axis equal
axis([min(rg) max(rg) min(zg) max(zg)])

