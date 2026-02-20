function [eqs, signals] = gsp_postprocessing(t, eqs, signals, tok, settings, L)
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

for i = 1:length(eqs)
  eqs{i}.t = t(i);
end

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
    strike_info = sparc_strikept_finder(eqs{i}, tok, 0);   
    eqs{i} = copyfields(eqs{i}, strike_info, fieldnames(strike_info), 0);
    for sig = sigs
      signals.(sig{:}).Data(i,:) = strike_info.(sig{:});
    end
  end
end

% Compute a few more plasma parameters for SPARC
if settings.calc_post_prc_ext
    % Compute the parallel connection length Lz [m]
    psinvals = 1.0001;  % corresponds to ~1 labmda_q = 0.3mm for 8.7MA equilbrium
    plotit = 0;
    Lz = calc_Lz_profile(eqs, tok, psinvals, plotit, 0, t);

    for i = 1:length(eqs)
        eqs{i}.Lz = Lz.Data(i); 

        % Compute the 2D toroidal flux vs (R, Z)
        % FtPQ is the total toroidal flux contained in flux surface
        if sum(eqs{i}.FtPQ)>0.0
            FtPx = interp1(eqs{i}.F0+(eqs{i}.F1-eqs{i}.F0)*(L.pQ).^2,eqs{i}.FtPQ,eqs{i}.Fx,'linear',0);
            FtPx(L.lxy) = FtPx(L.lxy).*(eqs{i}.Opy(:)>0);
            FtPx(~L.lxy) = 0;
            eqs{i}.FtPx = FtPx;
        else
            eqs{i}.FtPx = 0.0*L.lxy;
        end
    end
end

% Append boundary-related post-processing quantities to LY
if settings.boundary_post_run
    % Inputs: (L, LY, dr_OMP, direction, npts)
    % dr_OMP: radial distance from LCFS to field line start point 
    % on Outboard Mid Plane (OMP) (can be a vector, default 0.001 [m]).
    % direction: direction for field line tracing.
    % +1 defalut, means following the B direction
    % npts: number of points to down-interpolate field line trajectory to (default: 101)
    for i = 1:length(eqs)
        LY_bdy = meqpostboundary(L,eqs{i},settings.dr_OMP, settings.direction, settings.npts);
        eqs{i} = LY_bdy;
    end    
end

eqs = eqs2LY(eqs);

end

%% local functions

function LYs = eqs2LY(eqs)
  % eqs is the same as LY but has some additional fields
  % remove non-LY fields to match structure of default MEQ LY struct
  fds2remove = {
      'islimited',                 
      'outside_plasma_domain_mask',
      'pcurrt',                    
      'psiapp',                    
      'psibry',                    
      'psimag',                    
      'psipla',                    
      'psitouch',                  
      'psix',                      
      'psizr',                     
      'rbbbs',                     
      'rbdef',                     
      'rg',                        
      'rl',                        
      'rmaxis',                    
      'rtouch',                    
      'rx',                        
      'zbbbs',                     
      'zbdef',                     
      'zg',                        
      'zl',                        
      'zmaxis',                    
      'ztouch',                    
      'zx'                        
    };
  LYs = cell(length(eqs),1);
  for i = 1:length(eqs)
    LYs{i} = rmfields(eqs{i}, fds2remove);
  end
end

