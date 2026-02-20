function [Lp, Li, Le, li, mcIp, mvIp] = inductance(eq, tok_data_struct)
% *************************************************************************
% Description: 
%    Compute various inductance-related quantities
% 
% Inputs: 
%   eq - equilibrium object
%   tok_data_struct - Toksys geometry object
% 
% Outputs: 
%  Lp - total plasma inductance
%  Li - internal inductance
%  Le - external inductance
%  li - normalized internal inductance
%  mcIp - mutuals coupling between coils and total plasma
%  mvIp - mutuals coupling between vessel and total plasma
% *************************************************************************
eq.pcurrt = eq.pcurrt * eq.cpasma / sum(eq.pcurrt(:));

struct_to_ws(tok_data_struct); 
mpc = tok_data_struct.mpc;

i = eq.rbbbs==0 & eq.zbbbs==0;
eq.rbbbs(i) =[];
eq.zbbbs(i)=[];

[~,psi_r,psi_z] = bicubicHermite(rg, zg, eq.psizr, rgg, zgg);

Br = -1./(2*pi*rgg) .* psi_z;
Bz =  1./(2*pi*rgg) .* psi_r;

dr = mean(diff(rg));
dz = mean(diff(zg));

in = inpolygon(rgg, zgg, eq.rbbbs, eq.zbbbs);
Br(~in) = 0;
Bz(~in) = 0;

mu0 = pi*4e-7;

volumezr = 2*pi*rgg*dr*dz;

ip = sum(eq.pcurrt(:));


Wmag = 1 / (2*mu0) * sum(sum((Br.^2 + Bz.^2).*volumezr));  % magnetic field energy

Li = 2*Wmag / ip^2;

R0 = (min(eq.rbbbs) + max(eq.rbbbs))/2;

Vtot = sum(sum(volumezr(in)));

Bp2volavg = sum(sum((Br.^2 + Bz.^2).*volumezr)) / Vtot;
Cl =  sum(sqrt(diff(eq.rbbbs).^2 + diff(eq.zbbbs).^2));
Bp2bryavg = (mu0*ip/Cl)^2;

li = Bp2volavg / Bp2bryavg;

% several ways to estimate psizr_pla, depending on inputs
if isfield(eq, 'psizr_pla')  
  psizr_pla = eq.psizr_pla;
elseif isfield(eq, 'pcurrt')
  psizr_pla = mpp_x_vec(mpp, eq.pcurrt(:));
  psizr_pla = reshape(psizr_pla, nz, nr);
else
  psizr_app = reshape(mpc*eq.ic + mpv*eq.iv, nz, nr);
  psizr_pla = eq.psizr - psizr_app;
end
  
[rbbbs, zbbbs] = interparc(eq.rbbbs, eq.zbbbs, 100, 0, 0);
psibry_pla = bicubicHermite(rg, zg, psizr_pla, rbbbs, zbbbs);
Le = mean(psibry_pla) / ip;

% inductance between plasma and conductors
if isfield(eq, 'pcurrt')
  circ = sparc_circ(tok_data_struct);
  mcIp = circ.Pcc' * mpc' * eq.pcurrt(:) / ip;
  mvIp = circ.Pvv' * mpv' * eq.pcurrt(:) / ip;
else
  mcIp = nan;
  mvIp = nan;
end

Lp = Le + Li;