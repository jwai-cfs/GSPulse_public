function [eqs, signals] = gsp_postprocessing(t, eqs, signals, tok, settings)
% =========================================================================
% Description: 
%  Intended as a catch-all function for doing any post-processing on the
%  equilibrium, after the last GS-iteration. 
%
% Inputs: 
%  t - time base for equilibria
%  eqs - cell array of equilibria
%  tok - tokamak geometry object, see help _define_tok.m
%  L   - MEQ geometry object if available, else empty
%  settings - settings struct, see help _define_settings.m
%
% Outputs: 
%  eqs - cell array of equilibria, with boundaries traced and available in
%    (eqs{}.rbbbs, eqs{}.zbbbs)
%
% =========================================================================

% trace boundary
if settings.do_final_boundary_trace
  eqopts = struct;
  eqopts.plotit = 0;
  eqopts.warmstart = [];
  eqopts.return_bry = 1;
  for i = 1:length(eqs)
    if settings.verbose 
      fprintf('    tracing boundary %d/%d, t=%.3f\n', i, length(eqs), t(i));
    end
    [eq, ws] = get_eqinfo(tok.rg, tok.zg, eqs{i}.psizr, tok.rl, tok.zl, eqopts);
    eqopts.warmstart = ws;
    eqs{i} = copyfields(eqs{i}, eq, [], 1);  
  end
end


% compute strike points (only SPARC supported)
if settings.calc_strike_pts && strcmpi(settings.tokamak, 'sparc')
  addpath([getenv('GSROOT') '/tokamaks/sparc/analysis'])  
  sigs = {'rstrike', 'zstrike', 'rtouch', 'ztouch', 'is_limited', 'is_UN', 'is_LN'};
  for sig = sigs
    signals.(sig{:}).Time = t;
  end
  for i = 1:length(eqs)
    eqs{i}.strike_info = sparc_strikept_finder(eqs{i}, tok, 0);   
    for sig = sigs
      signals.(sig{:}).Data(i,:) = eqs{i}.strike_info.(sig{:});
    end
  end
end






