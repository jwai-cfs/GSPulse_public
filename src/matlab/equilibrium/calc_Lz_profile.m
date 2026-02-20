function Lz = calc_Lz_profile(soln_input, tok, psinvals, plotit, debug_plots, varargin)
% =========================================================================
% Compute connection length Lz profile for each equilibria in the gspulse 
% soln. Lz is computed at specific values of psiN (e.g. for psiN = 
% 1.001: 0.001: 1.01) which are specified by psinvals. 
%
% EXAMPLE:
% soln = load('soln103.mat').soln;
% settings.nvessmodes = 40;
% [~, tok] = define_tok(settings);
% psinvals = 1.002: 0.002: 1.01;
% plotit = 1;
% Lz = calc_Lz_profile(soln, tok, psinvals, plotit);
% =========================================================================

if ~exist('plotit', 'var'), plotit = 0; end
if ~exist('debug_plots', 'var'), debug_plots = 0; end

Lz = struct;
if isempty(varargin)
    neq = length(soln_input.eqs);
    Lz.Time = soln_input.t;
else
    neq = length(soln_input);
    Lz.Time = varargin{1};
end

npts = length(psinvals);
Lz.Data = nan(neq, npts);
rg = tok.rg;
zg = tok.zg;
rlim = tok.limdata(2,:);
zlim = tok.limdata(1,:);
rBt = 1.85 * 12.2;

for i = 1:neq

  % read equilibrium
  if isempty(varargin)
      eq = soln_input.eqs{i};
  else
      eq = soln_input{i};
  end
  psizr = eq.psizr;
  if ~(isempty(eq.FA) || isempty(eq.FB))
      psinzr = (eq.psizr - eq.FA) / (eq.FB - eq.FA);
      seg = [eq.rA eq.rA+1 eq.zA eq.zA];
    
      % find starting points along outer midplane segment
      opts = struct;
      opts.plotit = debug_plots;
      [rs, zs] = seg_strike_finder(rg, zg, psinzr, psinvals, seg, opts);
    
      % calc connection length for each
      doplot = 1;
      for j = 1:npts
        Lz.Data(i,j) = calc_Lz(rg, zg, -psizr, rs(j), zs(j), rBt, eq.rA, ...
          eq.zA, rlim, zlim, debug_plots);
      end
  end
end

% plotting
if plotit
  figure
  plot(Lz.Time, Lz.Data, 'linewidth', 1.5)
  title('Connection Length', 'fontsize', 16)
  xlabel('Time [s]', 'fontsize', 14)
  ylabel('Lz [m]', 'fontsize', 14)
  legend(strcat('\psi_N =', num2str(psinvals')), 'fontsize', 12, 'location', 'northwest')
end
end

