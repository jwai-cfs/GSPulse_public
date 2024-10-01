function [eqs, pcurrt] = gs_update_psipla(...
  kinterval, mpcsoln, pcurrt_prev, tok, plasma_params, settings)
% =========================================================================
% Description: 
% Given the coil and vessel currents and previous equilibrium solution,
% calls out to FBT to update the plasma current and flux distribution of
% the equilibrium.
% 
% Inputs: 
%  kinterval - which time interval is being solved for
%  mpcsoln - solution of the MPC quadratic program from mpc_update_psiapp.m
%  pcurrt_prev - waveform for the previous iteration evolution of plasma
%                current
%  tok - tokamak geometry object, see help _define_tok.m
%  plasma_params - plasma parameters object, see help _define_plasma_params.m
%  settings - settings object, see help _define_settings.m
% 
% Outputs: 
%  eqs - cell array, next iteration of equilibrium
%  pcurrt - waveform for the updated evolution of plasma current
%
% Additional info: 
%  A written description of this part of the algorithm is given in
%  "Appendix B Step 5: Grad-Shafranov Picard iteration" in the
%  gspulse_algorithm.pdf. Equation references in the code refer to this
%  document. 
%
% =========================================================================

psizr = mpcsoln.psizr.Data';
psiapp = mpcsoln.psiapp.Data';


% read inputs
timedat = settings.timedata.interval(kinterval);
N = timedat.N;
t = timedat.t;

% initialize
pcurrt.Time = t;
pcurrt.Data = pcurrt_prev.Data * nan;
eqs = cell(N,1);
eqopts = struct;
eqopts.plotit = 0;
eqopts.warmstart = [];
eqopts.return_bry = 0;
eqs0 = cell(N,1);

% define grid info
mu0 = pi*4e-7;
[rgg, ~] = meshgrid(tok.rg, tok.zg);
dr = mean(diff(tok.rg));
dz = mean(diff(tok.zg));
dA = dr*dz;

% Loop through all equilibria and update
for i = 1:N

  if settings.verbose    
    fprintf('    plasma picard %d/%d, t=%.3f\n', i, N, t(i));
  end

  % analyze the equilibrium from the flux distribution (find xpoints,
  % touchpoints, magnetic axis, etc)
  [eqs0{i}, ws] = get_eqinfo(tok.rg, tok.zg, psizr(:,i), tok.rl, tok.zl, eqopts);
  eqopts.warmstart = ws;

  % read the plasma parameters for plasma current, thermal energy, and
  % profile basis shapes
  eq = eqs0{i};
  ref = structts2struct(plasma_params, {'Ip','Wk','pprime','ttprime'}, t(i));

  psinzr = (eq.psizr - eq.psimag) / (eq.psibry - eq.psimag);
  out = eq.outside_plasma_domain_mask;
  psin = linspace(0,1,length(ref.pprime))';

  b = struct;

  % pressure basis - eqn B.41
  b.pres = cumtrapz(psin, ref.pprime)';
  b.pres = b.pres-b.pres(end);
  b.preszr = interp1([psin; 1000], [b.pres; 0], psinzr(:));
  b.preszr(out) = 0;

  % pprime basis - eqn B.36
  b.pprime = ref.pprime * 2*pi/(eq.psibry-eq.psimag);
  b.pprimezr = interp1([psin; 1000], [b.pprime'; 0], psinzr(:));
  b.pprimezr(out) = 0;

  % ttprime basis - eqn B.36
  b.ttprime = ref.ttprime * 2*pi/(eq.psibry-eq.psimag);  
  b.ttprimezr = interp1([psin; 1.001; 1000], [b.ttprime'; 0; 0], psinzr(:));
  b.ttprimezr(out) = 0;

 
  % set up the equations:  [Ip; Wk] = H * [cp; cf];
  % then this can be solved for [cp; cf] the P' and FF' coefficients
  H = zeros(2,2);
  R = rgg(:);
  H(1,:) = [R'*b.pprimezr  1./(mu0 * R') * b.ttprimezr] * dA;  % eqn B.39
  H(2,1) = 3*pi*R'*b.preszr*dA;                                % eqn B.40 
  c = H \ [ref.Ip; ref.Wk];

  % find the plasma current distribution from Grad-Shafranov eqn B.2
  jphi = [R.*b.pprimezr  1./(mu0 * R).*b.ttprimezr] * c;
  jphi = reshape(jphi, tok.nz, tok.nr);
  pcurrt_i = jphi * dA;

  % write data to equilbrium
  eq.pcurrt = pcurrt_i;
  eq.psiapp = reshape(psiapp(:,i), tok.nz, tok.nr);
  eq.psipla = tok.mpp * pcurrt_i(:);
  eq.psipla = reshape(eq.psipla, tok.nz, tok.nr);
  eq.psizr = eq.psiapp + eq.psipla;
  eq.pprime = b.pprime * c(1);
  eq.ttprime = b.ttprime * c(2);
  eq.psin = psin;
  
  eqs{i} = eq;
  pcurrt.Data(i,:) = pcurrt_i(:);
end
pcurrt.Data = smoothdata(pcurrt.Data, 1, 'gaussian', 5);
