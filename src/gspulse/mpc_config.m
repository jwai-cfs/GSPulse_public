function config =  mpc_config(kinterval, tok, shapes, settings, weights, targs)
% =========================================================================
% Description: 
%  builds dynamic model and some of the large MPC optimization matrices 
%  that only need to be computed once
% 
% Inputs: 
%  kinterval - index of which time interval to build the config for
%  tok           - tokamak geometry struct, see help _define_tok.m
%  shapes        - shapes struct, see help _define_shapes.m
%  settings      - settings struct, see help _define_settings.m
%  weights       - weights struct, see help _define_weights.m
%  targs         - targets struct, see help _define_targets.m
% 
% Outputs: 
%  config - struct with lots of optimization-related quantities defined
%
% Additional info: 
% see Appendix B, Step 4: optimize conductor evolution in 
% gspulse_algorithm.pdf for more details on the computations
%
% =========================================================================
% read parameters
interval_time_data = settings.timedata.interval(kinterval);
dt = interval_time_data.dt;
t = interval_time_data.t;
Nlook = interval_time_data.N;
wts = weights.wts;
dwts = weights.dwts;
d2wts = weights.d2wts;
fds2control = settings.fds2control;
nv = settings.nvessmodes;

%% build dynamics and output models

% build the dynamics model A,B matrices: eqns B.8, B.9, and B.10
M = [tok.mcc tok.mcv; tok.mcv' tok.mvv];
M = (M + M') /  2;
R = diag([tok.resc; tok.resv]);                   
Minv = inv(M);
A = -Minv * R;
B = Minv(:,1:tok.nc);

% build the output C matrices: eqn B.13
dpsizrdx = [tok.mpc tok.mpv];
Cmats = output_model(dpsizrdx, tok, shapes, t, settings.fds2control);

% read dims
[nx, nu] = size(B);
ny = size(Cmats{1},1);
nc = nx - nv;

%%
cv = define_data_indices(settings, targs, tok);
mpc_dims_errchck(settings, shapes, weights, targs, Cmats, cv);

%% Build the MPC prediction model

% discretize model: eqn B.11

[Ad,Bd] = zoh_discretize(A,B,dt);

% Prediction model used in MPC
nw = nx;
E = [];
F  = [];
Fw = [];
Apow  = eye(nx);
F_row = zeros(nx, Nlook*nu);
Fw_row = zeros(nx, Nlook*nw);

for i = 1:Nlook

  % eqn B.23
  F = [F; F_row];
  idx = (nu*i+1):(nu*(i+1));
  F_row = Ad * F_row;
  F_row(:,idx) = Bd;
  
  % eqn B.23
  Fw = [Fw; Fw_row];
  idx = (nw*i+1):(nw*(i+1));
  Fw_row = Ad * Fw_row;
  Fw_row(:,idx) = eye(nx);
  
  % eqn B.23
  E = [E; Apow];
  Apow = Ad * Apow;
  
end
Chat = blkdiag(Cmats{:}); % eqn B.17
M = -Chat * [E F];        % eqn B.26

%% Weighting matrices

% build u-transform matrices: eqns B.18 & B.19
wu   = structts2vec(wts, {'v'}, t-dt);
wdu  = structts2vec(dwts, {'v'}, t(2:end)-dt);
wd2u = structts2vec(d2wts, {'v'}, t(2:end-1)-dt);

Su = speye(Nlook*nu);

m = vec2mat([-1 1], Nlook-1, 1);
Sdu = kron(m, eye(nu));
Sdu = sparse(Sdu);

m = vec2mat([1 -2 1], Nlook-2, 1);
Sd2u = kron(m, eye(nu));
Sd2u = sparse(Sd2u);

% build e-transform matrices: eqn B.20
we   = structts2vec(wts, fds2control, t);
wde  = structts2vec(dwts, fds2control, t(2:end));
wd2e = structts2vec(d2wts, fds2control, t(2:end-1));

Se = speye(Nlook*ny); 

m = vec2mat([-1 1], Nlook-1, 1);
Sde = kron(m, eye(ny));
Sde = sparse(Sde);

m = vec2mat([1 -2 1], Nlook-2, 1);
Sd2e = kron(m, eye(ny));
Sd2e = sparse(Sd2e);


% convert nans and empties to zero, these occur if solving for less than 3
% equilibria due to the timing assumptions
[we, wde, wd2e, Se, Sde, Sd2e] = nan_or_empty_to_zero(we, wde, wd2e, Se, Sde, Sd2e);
[wu, wdu, wd2u, Su, Sdu, Sd2u] = nan_or_empty_to_zero(wu, wdu, wd2u, Su, Sdu, Sd2u);

% build quadprog inputs: eqn B.21
He = Se'*diag(we)*Se + Sde'*diag(wde)*Sde + Sd2e'*diag(wd2e)*Sd2e;
He = sparse(He);

Hu = Su'*diag(wu)*Su + Sdu'*diag(wdu)*Sdu + Sd2u'*diag(wd2u)*Sd2u;
Hu = sparse(Hu);
Hx1 = zeros(nx);


%% compute move-blocking (spline basis) matrices
% for use as in eqns B.42 and B.43
if settings.use_spline_basis
  n_basis_pts = round(length(t(3:end)) * settings.spline_basis_ratio);
  t_basis = linspace(t(3), t(end), n_basis_pts);  
  [~,bm] = spline_basis(t_basis, nan(n_basis_pts,1), t(3:end));  
  bm = blkdiag(eye(2), bm);
  BM = kron(bm, eye(nu));
  BM = blkdiag(eye(nx), BM);  
  pinvBM = pinv(BM);
  BM = sparse(BM);
  pinvBM = sparse(pinvBM);
else
  BM = 1;  % identity blocking matrix
  pinvBM = 1;
end


%% save config
config = variables2struct(BM,pinvBM,Chat,nx,nu,nv,nc,A,B,Ad,Bd, ...
  E,F,Fw,Cmats,Hu,Hx1,He,M); 
