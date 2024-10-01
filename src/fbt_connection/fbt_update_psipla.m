% Given the coil and vessel currents and previous equilibrium solution,
% calls out to FBT to update the plasma current distribution and
% equilibrium.
%
function [eqs, pcurrt] = fbt_update_psipla(kiter, kinterval, mpcsoln, ...
  eqs_prev, pcurrt_prev, tok, plasma_params, settings, L)


% read inputs
timedat = settings.timedata.interval(kinterval);
N = timedat.N;
t = timedat.t;

% check that limiter is closed (otherwise can miss touch-points)
[L.G.rl, L.G.zl] = close_curve(L.G.rl, L.G.zl);
L.G.nl = length(L.G.rl);

% initialize
pcurrt.Time = t;
pcurrt.Data = pcurrt_prev.Data * nan;
eqs = cell(N,1);
LYt = cell(N,1);

% only do post-processing on last iteration
dopost = kiter == settings.niter;

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

% plasma picard update
for i = istart:N

  if settings.verbose    
    if dopost
      str = '    plasma picard and postprocessing (interval %d/%d) (equil %d/%d), t=%.3f\n';
    else
      str = '  plasma picard (interval %d/%d) (equil %d/%d), t=%.3f\n';
    end
    fprintf(str, kinterval, settings.timedata.n_intervals, i, N, t(i));
  end
  
  Iy = reshape(pcurrt_prev.Data(i,:), tok.nz, tok.nr);
  Iy = Iy(2:end-1, 2:end-1);
  Ia = structts2vec(mpcsoln, {'ic'}, t(i));
  Iu = structts2vec(mpcsoln, {'iv'}, t(i));

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

    % assign profile basis funs as expected by FBT code
    bfp = variables2struct(gNg, IgNg, fPg, fTg);
    L.P.bfp = bfp;
    L.bfp = bfp;
    L.fPg = bfp.fPg;
    L.fTg = bfp.fTg;

  end

  % get plasma current distro via fbt
  [Iy, LYt, Fx] = fbt_Iy_gspulse(L, Iy, Ia, Iu, LX, dopost);
  if ~dopost
    Fx = meqFx(L,Iy,[Ia;Iu]);
  end
  Iy = blkdiag(0, Iy, 0);
  pcurrt.Data(i,:) = Iy(:);

  eq = LYt;
  eq.pcurrt = Iy;
  eq.psiapp = tok.mpc*Ia + tok.mpv*Iu; 
  eq.psipla = Fx(:) - eq.psiapp;
  eq.psizr = Fx; 
  eq.outside_plasma_domain_mask = abs(eq.pcurrt) < sqrt(eps);

  eqs{i} = eq;
end

