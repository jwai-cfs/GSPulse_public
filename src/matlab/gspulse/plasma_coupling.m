function w_hat = plasma_coupling(pcurrt, tok, A)
% =========================================================================
% Description: 
%  compute the flux at the coils induced from the plasma current motion
%  
% Inputs: 
%  dt - time step for the pcurrt
%  pcurrt - plasma current array with dimension [nz*nr x length(time)]
%  tok           - tokamak geometry struct, see help _define_tok.m
%
% Outputs: 
%  w - plasma current motion coupling term
%
%
% Additional info: 
%   This term is described by eqns A.6, A.7, A.8, A.9, A.10 in
%   gspulse_algorithm.pdf
%  
% =========================================================================
M = [tok.mcc tok.mcv; tok.mcv' tok.mvv];
M = (M + M') /  2;                 
Minv = inv(M);

dt = diff(pcurrt.Time);
pcurrtdot = sparse(diag(1./dt)) * diff(pcurrt.Data);
w_hat_continuous = -Minv * [tok.mpc tok.mpv]' * pcurrtdot';

% discretize
w_hat = nan(size(w_hat_continuous));
for i = 1:length(dt)
  [~, w_hat(:,i)] = zoh_discretize(A, w_hat_continuous(:,i), dt(i));
end
