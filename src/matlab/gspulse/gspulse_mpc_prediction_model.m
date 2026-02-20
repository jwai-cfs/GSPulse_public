function [E,F,Fw,M] = gspulse_mpc_prediction_model(A, B, Cmats, Nlook, dt)

assert(Nlook == length(dt)+1)
assert(Nlook == length(Cmats))

% Initialize
[nx, nu] = size(B);
nw = nx;
ny = size(Cmats{1}, 1);
Bw = eye(nx);
E_row  = eye(nx);              % eqn A.23
F_row = zeros(nx, Nlook*nu);   % eqn A.23
Fw_row = zeros(nx, Nlook*nw);  % eqn A.23
M_row = -Cmats{1} * [E_row F_row]; % eqn A.26

E = zeros(Nlook*nx, nx);
F = zeros(Nlook*nx,Nlook*nu);
Fw = zeros(Nlook*nx, Nlook*nw);
M = zeros(Nlook*ny, nx + Nlook*nu);

% First timestep of prediction model
E(1:nx,:) = E_row;
F(1:nx,:) = F_row;
Fw(1:nx,:) = Fw_row;
M(1:ny,:) = M_row;

% Timesteps 2 to N of prediction model
for it = 2:Nlook
  
  [Ad,Bd] = zoh_discretize(A,B,dt(it-1));

  % eqn A.23
  i = (it-1)*nx + (1:nx);
  j = (nu*(it-1)+1):(nu*it);
  F_row = Ad*F_row;
  F_row(:,j) = Bd;
  F(i,:) = F_row;

  % eqn A.23
  j = (nw*(it-1)+1):(nw*it);
  Fw_row = Ad * Fw_row;
  Fw_row(:,j) = eye(nx);
  Fw(i,:) = Fw_row;
 
  % eqn A.23
  E_row = Ad * E_row;
  E(i,:) = E_row;
  
  % eqn A.26
  i = (it-1)*ny + (1:ny);
  M(i,:) = -Cmats{it} * [E_row F_row];

end

% Note: Implementation for Fw is amended from journal article. Article uses
% w_hat=[w0,w1,...,w_N-1] but we amend this to use w_hat=[w1,...,w_N-1] and 
% delete the corresponding columns of Fw. 
Fw(:,1:nw) = []; 
