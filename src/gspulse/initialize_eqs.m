function [eqs1N, t1N] = initialize_eqs(kinterval, settings, shapes, pcurrt, tok, init)
% =========================================================================
% Description: 
%  builds an initial estimate of the flux distribution of each equilibrium,
%  based on the crude estimate of plasma current from initialize_pcurrt.m
%  
% Inputs:  
%  kinterval - which time interval is being solved for
%  settings - settings struct, see help _define_settings.m
%  shapes - shapes struct, see help _define_shapes.m
%  pcurrt - waveform with "Time" and "Data" describing the plasma current
%         evolution estimate
%  tok - tokamak geometry struct, see help _define_tok.m
%  init - init struct, see help _define_init.m
% 
% Outputs: 
%  eqs1N - cell array with equliibria from t=1 to t=N
%  t1N   - timebase for the equilibria
% =========================================================================

% read inputs
timedat = settings.timedata.interval(kinterval);
t1N = timedat.t1N;

i = ~isnan(init.x1);
mpx = [tok.mpc tok.mpv];
psiapp1 = mpx(:,i) * init.x1(i);


psipla1N = tok.mpp * pcurrt.Data';
psizr1N = psipla1N + psiapp1; 
[rgg, zgg] = meshgrid(tok.rg, tok.zg);

% build dummy equilibria
eqs1N = {};
for i = 1:length(t1N)  
  ref = structts2struct(shapes, fieldnames(shapes), t1N(i));
  out = ~inpolygon(rgg, zgg, ref.rb, ref.zb);

  eq.psipla = psipla1N(:,i);
  eq.psizr = reshape(psizr1N(:,i), tok.nz, tok.nr);
  eq.outside_plasma_domain_mask = out;
  eq.psibry = bicubicHermite(tok.rg, tok.zg, eq.psizr, ref.rbdef, ref.zbdef);
  eqs1N{i} = eq;
end
