function [mpcsoln, wsout] = mpc_update_psiapp(kinterval, pcurrt, config, ...
  tok, shapes, init, settings, targs, model_offset, ws, eqs)
% =========================================================================
% Description: 
%  sets up and solves an MPC-like quadratic program for the dynamic
%  equilibrium optimization problem
%
% Inputs: 
%  kinterval - index of which time interval to build the config for
%  pcurrt    - waveform describing plasma current evolution
%  config    - mpc config object, output of mpc_config.m
%  tok           - tokamak geometry struct, see help _define_tok.m
%  shapes        - shapes struct, see help _define_shapes.m
%  init          - init struct, see help _define_init.m
%  settings      - settings struct, see help _define_settings.m
%  targs         - targets struct, see help _define_targets.m
%  model_offset  - model_offset struct, see help _define_model_offset.m
%  ws            - warmstart object for matlab's quadprog. Used to accelerate
%                  the solve after the first iteration
% 
% Outputs: 
%  mpcsoln - struct of timeseries containing waveforms for the voltages,
%            currents, shaping errors, etc
%  wsout   - warmstart object from matlab's quadprog.m
%
% Additional info: 
%  see Appendix B, Step 4: optimize conductor evolution in 
%  gspulse_algorithm.pdf for more details on the computations
%
% =========================================================================

% read parameters
timedat = settings.timedata.interval(kinterval);
Nlook = timedat.N;
t = timedat.t;
dt = timedat.dt;
fds2control = settings.fds2control;
nu = config.nu;
nx = config.nx;
E = config.E;
F = config.F;
Fw = config.Fw;
Chat = config.Chat;
Hu = config.Hu;
He = config.He;
Hx1 = config.Hx1;
M = config.M;
BM = config.BM;
pinvBM = config.pinvBM;
cv = define_data_indices(settings, targs, tok);

% read plasma flux
psipla = zeros(tok.nz*tok.nr, Nlook);
for i = 1:Nlook
  psipla(:,i) = eqs{i}.psipla(:);
end

% inject in the modelling offset
if settings.inject_model_offset
   psi_model_offset = structts2struct(model_offset, {'psizr'}, t).psizr';
else
  psi_model_offset = 0;
end

% measure y
ic = zeros(tok.nc, Nlook);
yks = measure_ys(psipla - psi_model_offset, ic, fds2control, shapes, tok, t);
ykhat = vertcat(yks{:});

% target y and error 
rhat = structts2vec(targs, fds2control, t);
dytarghat = rhat - ykhat;
ny = length(yks{1});

% plasma-motion induced coupling term: see eqns B.8, B.9, B.10, B.11
if isfield(init, 'pcurrt0')
  pcurrtext = [init.pcurrt0(:) pcurrt.Data'];
else
  pcurrtext = [zeros(tok.nz*tok.nr,1) pcurrt.Data'];
end
w = plasma_coupling(dt, pcurrtext, tok);

[~,wd] = zoh_discretize(config.A, w, dt);  
wd = wd(:);


% effect of plasma current motion on the QP optimization, eqn B.26
Fwwd = Fw * wd;
if length(t) == 1, Fwwd = zeros(size(Fwwd)); end
d = dytarghat - Chat * Fwwd;

% modify the prediction model when using multiple intervals, eqn B.45
if kinterval > 1  
  i = 1:(ny*2);  
  d(i) = [init.e1; init.e2];
  M(i,:) = 0;
end

% quadprog inputs: eqn B.28
f = M' * He * d;
H = blkdiag(Hx1, Hu) + M' * He * M;

% inequality constraints on the voltage (plus constraint on x1) 
if settings.enforce_voltage_limits
  ub_uhat = repmat(settings.vmax, Nlook, 1);
  lb_uhat = repmat(settings.vmin, Nlook, 1);  
  ub = [inf(nx,1); ub_uhat];
  lb = [-inf(nx,1); lb_uhat];  
else
  ub = inf(Nlook*nu + nx, 1);
  lb = -inf(Nlook*nu + nx, 1);
end

% current limits
if settings.enforce_current_limits  
  ymin = -inf(ny,1);
  ymax = inf(ny,1);
  ymin(cv.iy.ic) = settings.ic_min;
  ymax(cv.iy.ic) = settings.ic_max;
  yminhat = repmat(ymin, Nlook, 1);
  ymaxhat = repmat(ymax, Nlook, 1);  
  
  if kinterval > 1   
    % except for the first stage, constraints on the first 2 timesteps are
    % enforced separately (implicitly via e1/e2/x1/u0/u1), so we need to
    % remove the current constraint from the first 2 timesteps here or 
    % else the QP can be overconstrained and not pass numerical checks. 
    yminhat(1:2*ny) = -inf;
    ymaxhat(1:2*ny) = inf; 
  end

  Aineq = [-M; M];  % eqn B.31
  bineq = [ymaxhat + d - rhat; -yminhat - d + rhat];

  i = isinf(bineq);
  Aineq(i,:) = [];
  bineq(i,:) = [];
else
  Aineq = [];
  bineq = [];
end


% spline basis compression (move blocking): eqns B.42 and B.43
f = BM' * f;
if ~isempty(Aineq), Aineq = Aineq * BM; end
if ~isempty(ub), ub = pinvBM * ub; end
if ~isempty(lb), lb = pinvBM * lb; end
H = BM' * H * BM;
H = (H+H')/2;


% equality constraints on x1, u0, and u1
% the constraints are enforced via construction here, instead of via the 
% (Aeq, beq) in quadrog due to numerical issues with that approach. 
if length(init.x1)==1, init.x1 = init.x1 * ones(nx,1); end
if length(init.u0)==1, init.u0 = init.u0 * ones(nu,1); end
if length(init.u1)==1, init.u1 = init.u1 * ones(nu,1); end

% standard case
if length(t) > 1  
  iconst = 1:(nx+2*nu);   
else
  % if only solving for 1 time, then u1 does not enter the solution, and
  % the optimizer only solves for phat = [x1; u0]. 
  init.u1 = []; 
  iconst = 1:(nx+nu); 
end   

% initialize the optimizer solution phat: eqn B.26
np = Nlook*nu + nx;  % number of optimization variables
phatb = pinvBM * nan(np, 1);
phatb(iconst) = [init.x1; init.u0; init.u1];  

% separate the solution into free and locked variables, depending on which
% initial conditions were specified as nan
ifree = isnan(phatb);     % index of optimization vars that are free
ilocked = ~isnan(phatb);  % index of optimization vars that are locked (equality constraint)

Hf = H(ifree, ifree);
ff = f(ifree) + H(ifree, ilocked) * phatb(ilocked);
Aineqf = Aineq(:, ifree);
bineqf = bineq - Aineq(:,ilocked) * phatb(ilocked);
ubf = ub(ifree);
lbf = lb(ifree);


% warmstart the quadprog  
recalculate_ws = isempty(ws) ||  (length(ws.X) ~= length(ff)) || ...
                 any(ws.X - sqrt(eps) > ubf) || any(ws.X + sqrt(eps) < lbf);
if recalculate_ws
  x0 = ff*0;
  x0(x0 > ubf) = ubf( x0 >ubf);
  x0(x0 < lbf) = lbf( x0 <lbf);
  if isequal(settings.qpsolver,'quadprog')
    qpopts = optimoptions('quadprog', 'Algorithm','active-set','Display','off');
    ws = optimwarmstart(x0, qpopts);
  end
end

% solve quadprog
if settings.verbose
  fprintf('  solving coil current optimization (interval %d/%d) ...\n', ...
    kinterval, settings.timedata.n_intervals)
end

switch settings.qpsolver
  case 'quadprog'
    [wsout, ~, exitflag] = quadprog(Hf, ff, Aineqf, bineqf, [], [], lbf, ubf, ws);
  case 'ipm'
    ipmopts = {'debug',settings.verbose>0};
    % warn if we're ignoring bounds
    if numel(bineqf)>0, warning('ignoring inequality constraints when using ipm'); end
    hasbounds = (numel(lbf)>0 && any(~isinf(lbf))) || ...
      (numel(ubf) && any(~isinf(ubf)));
    if hasbounds, warning('ignoring lower/upper bounds when using ipm'); end
    % remove bounds
    Aineqf = []; bineqf = []; lbf = []; ubf = [];
    
    [wsout.X,~,exitflag,~,~] = ipmwrapper(Hf,ff,Aineqf,bineqf,[],[],lbf,ubf,x0, ...
      ipmopts{:});
  otherwise
    error('not supported qpsolver %s',settings.qpsolver)
end
phatb(ifree) = wsout.X;

if exitflag ~= 1
  str = sprintf('quadprog did not converge with exitflag %d\n', exitflag);
  warning(str);
end

% de-transform the move blocking
phat = BM * phatb;  % eqn B.42
uhat = phat(nx+1:end);

% extract predictions
xhat = [E F] * phat + Fwwd;  % eqn B.24
ehat = M * phat + d;         % eqn B.25
yhat = rhat - ehat;          % eqn B.25 

% save predictions with structts data format
ydat = vec2structts(yhat, fieldnames(cv.iy), cv.iy, t);
edat = vec2structts(ehat, fieldnames(cv.ie), cv.ie, t);
udat = vec2structts(uhat, fieldnames(cv.iu), cv.iu, t-dt);
xdat = vec2structts(xhat, fieldnames(cv.ix), cv.ix, t);
xdat.iv = xdat.ivb;

% flux reconstructions
psiapp = [tok.mpc tok.mpv] * [xdat.ic.Data'; xdat.iv.Data'];
psiapp = psiapp - psi_model_offset;
psizr = psiapp + psipla;
ydat.psizr.Time = t;
ydat.psiapp.Time = t;
ydat.psipla.Time = t;
ydat.psizr.Data = psizr';
ydat.psiapp.Data = psiapp';
ydat.psipla.Data = psipla';

% consolidate 
mpcsoln = ydat;
mpcsoln = copyfields(mpcsoln, edat, [], 0);
mpcsoln = copyfields(mpcsoln, xdat, [], 0);
mpcsoln = copyfields(mpcsoln, udat, [], 0);
