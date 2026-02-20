function [Tuv, Tvu] = vvbal(mcc, mcv, mvv, resv, nu)
% =========================================================================
% Description:
% produces balancing transformation matrices Tuv and Tvu such that: 
%
% Iu = Tuv * Iv   (exact)
% Iv ~= Tvu*Iu    (approximate)
% 
% where Iu represents the mode currents. 
% 
% The transformation is produced from a balanced realization of the 
% coil-vessel dynamical system:
%
% Ivdot = -inv(mvv) * rv * Iv - Mvc * inv(mcc) * Icdot
%
% Inputs: mcc - coil to coil mutuals, mcv - coil to vessel mutuals, mvv - 
% vessel mutuals, resv - vessel resistances, nu - number of vessel modes to
% retain. 
%
%  Outputs: 
%  (Tuv, Tvu) - balancing transform matrices as in the description
% =========================================================================
A = -inv(mvv) * diag(resv);
B = - mcv' * inv(mcc);
C = eye(size(mvv));

P = ss(A,B,C,0);
[~,~,Tuv,~] = balreal(P);
Tvu = inv(Tuv);
Tuv = Tuv(1:nu,:);
Tvu = Tvu(:,1:nu);