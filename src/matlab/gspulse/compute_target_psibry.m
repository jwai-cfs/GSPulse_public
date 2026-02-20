function target_psibry = compute_target_psibry(eqs0N, t0N, shapes, plasma_params, tok)
% =========================================================================
% Description: 
%  Integrate the Ejima equation (eqn A.4 in gspulse_algorithm.pdf) to find
%  the boundary flux evolution
%  
% Inputs: 
%  eqs0N - cellarray with equilibria from t=0 to t=N
%  t0N   - time vector for the equilibria
%  shapes - shapes object, see help _define_shapes.m
%  plasma_params - plasma parameters object, see help _define_plasma_params.m
%  tok - tokamak geometry object, see help _define_tok.m
% 
% Outputs: 
%  target_psibry - waveform with fields "Time" and "Data" describing the  
%                  boundary flux evolution
%
% Additional info: 
%  Note that Wmag = 1/2 Li * Ip^2. Wmag is used here, but the term in 
%    1/2 Li * Ip^2 is used in eqn A.4. 
% =========================================================================

assert(isuniform(t0N), 'Ejima equation with nonuniform time not yet supported.')

% get time info
dt = mean(diff(t0N));
t = t0N(2:end);
t0 = t0N(1);

% get initial condition for psibry
ref = structts2struct(shapes, {'cp_r', 'cp_z'}, t0);
psibry0 = mean(bicubicHermite(tok.rg, tok.zg, eqs0N{1}.psizr, ref.cp_r, ref.cp_z));

% measure stored magnetic energy
Wmag = [];
for i = 1:length(eqs0N)
  Wmag(i,1) = get_Wmag(eqs0N{i}.psizr, tok, eqs0N{i}.outside_plasma_domain_mask);
end

% Ejima equation
Wmagdot = diff(Wmag) / dt;
Rp = structts2vec(plasma_params, {'Rp'}, t-dt/2);
Ip = structts2vec(plasma_params, {'Ip'}, t-dt/2);
psibrydot = -Rp.*Ip - 1./Ip .* Wmagdot;

% integrate
target_psibry.Time = t0N; 
target_psibry.Data = psibry0 + cumsum([0; dt*psibrydot]);

end

% ==================
% Function: get_Wmag
% ==================
function Wmag = get_Wmag(psizr, tok, out)
  psizr = reshape(psizr, tok.nz, tok.nr);
  mu0 = pi*4e-7;
  dr = mean(diff(tok.rg));
  dz = mean(diff(tok.zg));
  volumezr = 2*pi*tok.rgg*dr*dz;
  [psi_r, psi_z] = gradient(psizr, dr, dz);
  Br = -1./(2*pi*tok.rgg) .* psi_z;
  Bz =  1./(2*pi*tok.rgg) .* psi_r;
  Br(out) = 0;
  Bz(out) = 0;
  Wmag = 1 / (2*mu0) * sum(sum((Br.^2 + Bz.^2).*volumezr));  % magnetic field energy
end
