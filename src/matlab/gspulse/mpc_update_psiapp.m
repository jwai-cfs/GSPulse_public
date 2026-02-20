function [mpcsoln, wsout] = mpc_update_psiapp(kinterval, pcurrt, config, ...
  tok, init, settings, optimization_signals, ws, eqs, psibry_target_override)
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
% ----------------
timedat = settings.timedata.interval(kinterval);
Nlook = timedat.N;
t = timedat.t1N;
dt = timedat.dt1N;
nu = config.nu;
nx = config.nx;
E = config.E;
F = config.F;
Fw = config.Fw;
Chat = config.Chat;
M = config.M;
Sm = config.Sm;
Hn = config.Hn; 
Aineqn = config.Aineqn;
Aineqs = config.Aineqs;
cv = define_data_indices(optimization_signals, tok);
tol = settings.qpsolver_tol;
[~,~,~,targs] = assemble_weights_and_targets(optimization_signals);

if ~isempty(psibry_target_override)
  targs.psibry = psibry_target_override;
end


% update plasma measurements
% ---------------------------

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

% measure y (from plasma contribution only)
ic = zeros(tok.nc, Nlook); % zeros b/c plasma contribution to coil current is none
yks = measure_ys(psipla - psi_model_offset, ic, tok, t, optimization_signals);
ykhat = vertcat(yks{:});

% target y and error 
rhat = structts2vec(targs, cv.y_names, t);
dytarghat = rhat - ykhat;
ny = length(yks{1});

% plasma-motion induced coupling term: see eqns A.8, A.9, A.10, A.11
% Note: Implementation is amended from journal article. Article defines 
% w_hat = [w0, w1, ..., w_N-1] but note that w0 does not affect the
% optimization because the corresponding entries of Fw are zero. This 
% implementation uses w_hat = [w1, w2, ..., w_N-1] and adjusts Fw 
% correspondingly, so that there is no need to compute w0. 
w_hat = plasma_coupling(pcurrt, tok, config.A);

% effect of plasma current motion on the QP optimization, eqn A.26
Fw_w_hat = Fw * w_hat(:);
if isscalar(t), Fw_w_hat = zeros(size(Fw_w_hat)); end

% update QP inputs
% -----------------
d = dytarghat - Chat * Fw_w_hat;  % eqn A.29

% modify the prediction model when using multiple intervals
if kinterval > 1  
  i = 1:(ny*2);  
  d(i) = [init.e1; init.e2]; % eqn A.45
end

% quadprog inputs
f = M' * config.He' * d;  %  eqn A.31

% update beq (RHS of equality constraints), see get_equality_constraint.m
i = config.equality_constraint_data.idx_current_locks_t2N;
b = config.equality_constraint_data.b;
b{1} = -d(i);  % equality constraint version of eqn A.34
if kinterval > 1 
  b{2} = init.x1;
  if Nlook == 1
    b{3} = init.u0;
  else
    b{3} = [init.u0; init.u1];
  end
end
beq = vertcat(b{:});

% update bineq (RHS of inequality constraints), see get_inequality_constraint.m
bineq = config.inequality_constraint_data.bineq;
bineq(1:2*ny*Nlook) = [ config.inequality_constraint_data.ymaxhat + d - rhat; 
                       -config.inequality_constraint_data.yminhat - d + rhat];

% remove infs
i = isinf(bineq);
bineq(i) = [];
Aineqn(i,:) = [];
Aineqs(i,:) = [];

% transformations
% -------------------

% spline transformation
fs = Sm' * f;

% equality-constraint-direct-substitution transformation
ceqs = config.pinvAeqs * beq; 
fn = config.nullAeqs' * (config.Hs' * ceqs + fs);
bineqn = bineq - Aineqs*ceqs;


% solve the QP and perform back-transformations
% ---------------------------------------------

% solve QP
wsout = solve_qp(Hn, fn, Aineqn, bineqn, settings.qpsolver, settings.verbose, tol, ws);

% equality-constraint-direct-substitution reverse transformation
chat = ceqs + config.nullAeqs * wsout.X; 

% spline reverse transformation (eqn A.45)
phat = Sm * chat;


% unravel prediction model
% --------------------------
uhat = phat(nx+1:end);

% extract predictions
xhat = [E F] * phat + Fw_w_hat;  % eqn A.24
ehat = M * phat + d;         % eqn A.25
yhat = rhat - ehat;          % eqn A.25 

% save predictions with structts data format
ydat = vec2structts(yhat, fieldnames(cv.iy), cv.iy, t);
edat = vec2structts(ehat, fieldnames(cv.ie), cv.ie, t);
udat = vec2structts(uhat, fieldnames(cv.iu), cv.iu, timedat.t0Nm1);
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
