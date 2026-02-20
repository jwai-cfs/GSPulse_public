function config =  mpc_config(kinterval, init, settings, tok, optimization_signals)
% =========================================================================
% Description: 
%  builds dynamic model and some of the large MPC optimization matrices 
%  that only need to be computed once
% 
% Inputs: 
%
%
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
dt = interval_time_data.dt1N;
t = interval_time_data.t1N;
Nlook = interval_time_data.N;
nv = length(tok.resv);
hpv = settings.verbose > 2; % hyper verbose

%% build dynamics and output models

% build the dynamics model 
% A,B matrices: eqns A.8, A.9, and A.10
Mut = [tok.mcc tok.mcv; tok.mcv' tok.mvv];
Mut = (Mut + Mut') /  2;
R = diag([tok.resc; tok.resv]);                   
Minv = inv(Mut);
A = -Minv * R;
B = Minv(:,1:tok.nc);

% build the output C matrices: eqn A.13
dpsizrdx = [tok.mpc tok.mpv]; 
vprintf(hpv,'building Cmat... '); tti=tic;
Cmats = output_model(dpsizrdx, tok, t, optimization_signals);
vprintf(hpv,' %f\n',toc(tti))

[nx,nu] = size(B);
ny = size(Cmats{1},1);
Chat = sparse(blkdiag(Cmats{:}));

% build the MPC prediction model: eqn A.23, A.26
vprintf(hpv,'building MPC matrices... '); tti=tic;
[E,F,Fw,M] = gspulse_mpc_prediction_model(A, B, Cmats, Nlook, dt);
vprintf(hpv,' %f\n',toc(tti));

% modify the prediction model when using multiple intervals, eqn A.45
if kinterval > 1  
  M(1:(2*ny),:) = 0;
end

% Weighting matrices
cv = define_data_indices(optimization_signals, tok);
[Hp,He] = build_H_weights(cv, Nlook, nx, nu, ny, interval_time_data, optimization_signals);
H = Hp + M' * He * M;

% Spline transform
Sm = get_spline_transform(t, settings, nu, nx);

equality_constraint_data = get_equality_constraint(...
  cv, optimization_signals,interval_time_data,init,nu,nx,ny,Nlook,M,kinterval);
Aeq = equality_constraint_data.Aeq;

inequality_constraint_data  = get_inequality_constraint(...
  settings, ny, nu, nx, Nlook, M, cv, kinterval);
Aineq = inequality_constraint_data.Aineq;

% spline transform for the qp problem (H,f,Aineq,bineq,Aeq,beq)
% --------------------------------------------------------------
Hs = Sm' * H * Sm;
% fs = Sm' * f;   % f is updated at each iteration, not available now
Aeqs = Aeq * Sm;   
Aineqs = Aineq * Sm;
% beqs: do nothing (beqs = beq)
% bineqs: do nothing (bineqs = bineq)

% equality-constraint-direct-substitution transforms of the qp
% see: mathworks.com/matlabcentral/fileexchange/183056-quadratic-program-direct-substitution
% ---------------------------------------------------------------
[nullAeqs, pinvAeqs] = null_and_pinv(Aeqs);
Hn = nullAeqs' * Hs * nullAeqs;
% fn: transform is fn = NA' * (Hs' * pinvAeqs * beq + fs);
%     cannot solve now, since beq changes per iteration
Aineqn = Aineqs * nullAeqs;
% bineqn: do nothing (bineqn = bineqs)
% Aeqn = [];  % no longer used
% beqn = [];  % no longer used


% Plot relative scaling of weights
if settings.plotlevel > 1
  scale_e = He * blkdiag(Cmats{:}) *  vecnorm(F,2,2);
  scale_u = sum(Hu,2);
  e = vec2structts(scale_e, cv.y_names, cv.iy, t);
  u = vec2structts(scale_u, cv.u_name, cv.iu, t);
  s = merge_structs(e, u);
  plot_structts(s, fieldnames(s), [], 3);
  sgtitle('Relative scaling of weights')
end

%% save config
config = variables2struct(nx,nu,nv,A,B,E,F,Fw,Cmats,M,Sm,Chat,He,...
  H, Hs, Hn, Aeqs, Aineqs, Aineqn, nullAeqs, pinvAeqs, ...
  equality_constraint_data, inequality_constraint_data);

end


%% Local functions
function vprintf(doprint,varargin)

if doprint
  fprintf(varargin{:})
end
end
