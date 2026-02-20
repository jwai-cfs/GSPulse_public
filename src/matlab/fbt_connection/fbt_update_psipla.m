% Given the coil and vessel currents and previous equilibrium solution,
% calls out to FBT to update the plasma current distribution and
% equilibrium.
%
function [eqs, pcurrt] = fbt_update_psipla(kiter, kinterval, mpcsoln, ...
  eqs_prev, pcurrt_prev, tok, plasma_params, settings, L, converged)

% read inputs
timedat = settings.timedata.interval(kinterval);
N = timedat.N;
t = timedat.t1N(:);

% initialize
pcurrt.Time = t;
pcurrt.Data = pcurrt_prev.Data * nan;
eqs = cell(N,1);
LYt = cell(N,1);

% only do post-processing on last iteration
dopost = kiter == settings.niter || converged;

% When stitching together solutions across multiple intervals, the
% equilibria are overlapped by 2 equilibria. That is eqs{1} and eqs{2} are
% derived from the previous stage. Need to write-in the data for these
% equilibria. If they are updated in the Picard loop, then they suffer
% from artificial drift in the vertical position, similar to static
% "forward" equilibrium solvers, since the currents for these 2 equilibria
% are locked by the MPC optimization and cannot do "feedback" compensation. 
if kinterval > 1
  eqs{1} = eqs_prev{1};
  eqs{2} = eqs_prev{2};
  pcurrt.Data(:,1) = pcurrt_prev.Data(:,1);
  pcurrt.Data(:,2) = pcurrt_prev.Data(:,2);
  istart = 3;
else
  istart = 1;
end

% Apply time-based smoothing to pcurrt - can help with convergence
pcurrt_prev.Data = nonuniform_movmean(pcurrt_prev.Time, ...
  pcurrt_prev.Data, settings.Iy_smooth_dt);

% find indices where plasma current density distribution will be updated 
% per Grad-Shafranov, defined by where Ip is greater than a threshold value
Ip = structts2vec(plasma_params, {'Ip'}, t);
idx_gs = find(Ip >= settings.fbt_Ip_threshold);
idx_gs(idx_gs < istart) = [];

% find indices where we will use simplified interpolation (not
% Grad-Shafranov) to update the plasma current. 
idx_non_gs = setdiff(istart:N, idx_gs);

% Loop through Grad-Shafranov updates via fbt
% -------------------------------------------
for k = 1:length(idx_gs)
  i = idx_gs(k);

  if settings.verbose > 1 
    fprintf('    plasma Grad-Shafranov picard update (equil %02d/%02d), t=%.3f\n', k, length(idx_gs), t(i));
  end

  % read inputs
  [L, Iy, Ia, Iu, LX] = parse_fbt_inputs(mpcsoln, pcurrt_prev, tok, ...
    plasma_params, L, t, i);
  

  % call picard update from fbt
  Iy_prev = Iy;
  [Iy, LYt, Fx] = fbt_Iy_gspulse(L, Iy, Ia, Iu, LX, dopost);

  % apply relaxation factor - can help with convergence
  Iy = settings.Iy_relax_factor*Iy + (1-settings.Iy_relax_factor)*Iy_prev;

  if ~dopost
    Fx = meqFx(L,Iy,[Ia;Iu]);
  end

  % write to pcurrt and eqs
  Iy = blkdiag(0, Iy, 0);
  pcurrt.Data(i,:) = Iy(:);  
  eqs{i} = write_to_eq(LYt, Iy, Fx, tok, Ia, Iu);
  
end

% Loop through non-Grad-Shafranov updates 
% ---------------------------------------
for k = 1:length(idx_non_gs)
  i = idx_non_gs(k);
  
  if settings.verbose > 1 
    fprintf('    non-GS picard update (equil %02d/%02d), t=%.3f\n', k, length(idx_non_gs), t(i));
  end
  
  % read inputs
  [L, ~, Ia, Iu, LX] = parse_fbt_inputs(mpcsoln, pcurrt_prev, tok, ...
    plasma_params, L, t, i);

  % get Iy from interpolating to nearest real GS equilibrium
  [~,k] = min(abs(idx_gs-i));
  idx_nearest = idx_gs(k);   
  Iy_nearest = eqs{idx_nearest}.pcurrt;
  Iy_zero = zeros(tok.nz, tok.nr); 
  alpha = Ip(i) / Ip(idx_nearest);
  Iy = alpha * Iy_nearest + (1-alpha) * Iy_zero;  

  % compute Fx
  Fx = meqFx(L,Iy(2:end-1,2:end-1),[Ia;Iu]);

  % spoof LYt 
  LYt = struct('shot',0,'t',0,'aq',[],'aW',[]);  
  if dopost
    ag = nan(3,1);    
    isaddl = 1;
    rBt = LX.rBt;
    
    [rA,zA,FA,dr2FA,dz2FA,drzFA,rX,zX,FX,dr2FX,dz2FX,drzFX,...
      rB,zB,FB,lB,lX,Opy,F0,F1] = meqpdom(Fx,LX.Ip,isaddl,L);

    LYt = meqpost(L,LYt,ag,...
      Fx,FA,FB,rA,zA,dr2FA,dz2FA,drzFA,rB,zB,lB,lX, ...
      rX,zX,FX,dr2FX,dz2FX,drzFX, ...
      rBt,Ia,Iu,Iy(2:end-1,2:end-1),Opy,F0,F1);
  end

  % write to pcurrt and eqs
  pcurrt.Data(i,:) = Iy(:);
  eqs{i} = write_to_eq(LYt, Iy, Fx, tok, Ia, Iu);

end



end

function [L, Iy, Ia, Iu, LX] = parse_fbt_inputs(...
    mpcsoln, pcurrt_prev, tok, plasma_params, L, t, i)

  % read currents
  Iy = reshape(pcurrt_prev.Data(i,:), tok.nz, tok.nr);
  Iy = Iy(2:end-1, 2:end-1);
  Ia = structts2vec(mpcsoln, {'ic'}, t(i));
  Iu = structts2vec(mpcsoln, {'iv'}, t(i));
  
  % read fbt plasma parameters
  LX = structts2struct(plasma_params, fieldnames(plasma_params), t(i));

  % assign basis funs
  if isequal(L.P.bfct, @bfabmex)
    L.P.bfct = @bfabmex;
    L.P.bfp  = [1 2];
    L.bfct = @bfabmex;
    L.bfp = [1 2];

  elseif isequal(L.P.bfct, @bf3imex)

    % read profile basis functions
    x = structts2struct(plasma_params, {'pprime', 'ttprime1', 'ttprime2'}, t(i));
    gNg = [x.pprime(:) x.ttprime1(:) x.ttprime2(:)];
    fPg = [1 0 0]';
    fTg = [0 1 1]';
    IgNg = bfprmex(gNg);

    % if provided explicitly, overwrite pressure profile computed by bfprmex
    if isfield(plasma_params, 'press')
      IgNg(:,1) = structts2vec(plasma_params, {'press'}, t(i));
    end

    % write profile basis functions as expected by fbt
    bfp = variables2struct(gNg, IgNg, fPg, fTg);
    L.P.bfp = bfp;
    L.bfp = bfp;
    L.fPg = bfp.fPg;
    L.fTg = bfp.fTg;

  end
end

function eq = write_to_eq(LYt, Iy, Fx, tok, Ia, Iu)
  eq = LYt;
  eq.pcurrt = Iy;
  eq.psiapp = tok.mpc*Ia + tok.mpv*Iu; 
  eq.psipla = Fx(:) - eq.psiapp;
  eq.psizr = Fx; 
  eq.outside_plasma_domain_mask = abs(eq.pcurrt) < sqrt(eps);  
end