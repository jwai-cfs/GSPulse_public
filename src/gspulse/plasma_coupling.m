function w = plasma_coupling(dt, pcurrt, tok)
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
%   This term is described by eqns B.6, B.7, B.8, B.9, B.10 in
%   gspulse_algorithm.pdf
%  
% =========================================================================
M = [tok.mcc tok.mcv; tok.mcv' tok.mvv];
M = (M + M') /  2;                 
Minv = inv(M);
pcurrtdot = diff(pcurrt, 1, 2) ./ dt;
w = -Minv * [tok.mpc tok.mpv]' * pcurrtdot;