function [Ad,Bd] = zoh_discretize(A,B,dt)
% Discretize state-space model matrices A,B
% 
% zero order hold: use trick that e^([A,B;0;0]*dt) = [Ad,Bd;0,I]
% Chi-Tsong Chen (1984). Linear System Theory and Design p125
% This has been verified against [Ad, Bd] = c2d(A,B,dt); 

nx = size(A,2); nu = size(B,2);
AB = [A,B]; AB = [AB;zeros(nu,nu+nx)];
ABd = expm(AB*dt);
Ad = ABd(1:nx,1:nx); Bd = ABd(1:nx,1+nx:end);
