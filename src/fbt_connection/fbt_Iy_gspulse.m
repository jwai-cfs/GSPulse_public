% gspulse version of fbt_Iy.m script from meq
%
% NOT FOR OPEN-SOURCE RELEASE

function [Iy, LYt, Fx] = fbt_Iy_gspulse(L, Iy, Ia, Iu, LX, dopost)

% configure inputs
shot = 0;
t = 0;
isaddl = 1;
ag = [0 0 0]';
Ip = LX.Ip;
rBt = LX.rBt;

%% Update flux distribution
Fx = meqFx(L,Iy,[Ia;Iu]);

%% Plasma domain and current distribution
[rA,zA,FA,dr2FA,dz2FA,drzFA,rX,zX,FX,dr2FX,dz2FX,drzFX,...
  rB,zB,FB,lB,lX,Opy,F0,F1,status,msg,id,dF0dFx,dF1dFx,ixI] = meqpdom(Fx,Ip,isaddl,L);


%% Internal plasma profile constraints
% Prerequisites
[Tyg,TpDg,ITpDg] = L.bfct(1,L.bfp,Fx,F0,F1,Opy,L.ry,L.iry);

% Fit using meqagcon
%  This amounts to a Newton iteration in ag space
[res, dresdF0,dresdF1,dresdag,dresdFx,dresdTpDg,dresdITpDg,...
          dresdrA, dresddr2FA, dresddz2FA, dresddrzFA, dresdCo] = ...
          meqagcon(L,LX,F0,F1,rA,dr2FA,dz2FA,drzFA,ag,Fx,Opy,TpDg,ITpDg);

ag = ag-dresdag\res;

%% plasma current from basis function coefficients
% Compute current density
Iy(:) = Tyg*ag;


%% Post-processing
aq = []; aW = [];  % no initial guess as time slices are independent
LYt = struct('shot',shot,'t',t,'aq',aq,'aW',aW);
if dopost
  LYt = meqpost(L,LYt,ag,...
  Fx,FA,FB,rA,zA,dr2FA,dz2FA,drzFA,rB,zB,lB,lX, ...
  rX,zX,FX,dr2FX,dz2FX,drzFX, ...
  rBt,Ia,Iu,Iy,Opy,F0,F1);
end