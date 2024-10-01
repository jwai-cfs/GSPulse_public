function pcurrt = initialize_pcurrt(tok, shapes, plasma_params, t)
% =========================================================================
% Description: 
%  Get a rough estimate of plasma current distribution for each shape, 
%  based on a simple parabolic distribution estimate of the current density
%  (eqn B.5)
%  
% Inputs:  
%  tok           - tokamak geometry struct, see help _define_tok.m
%  shapes        - shapes struct, see help _define_shapes.m
%  plasma_params - plasma parameters, see help _define_plasma_params.m
%  t - timebase on which to initialize
% 
% Outputs: 
%  pcurrt - waveform with "Time" and "Data" describing the plasma current
%  evolution estimate
%
% =========================================================================
opts.plotit = 0;

N = length(t);
nr = tok.nr;
nz = tok.nz;
pcurrt.Time = t;
pcurrt.Data = zeros(N, nz*nr);

for i = 1:N  
  ip = structts2vec(plasma_params, {'Ip'}, t(i));
  rb = structts2vec(shapes, {'rb'}, t(i));
  zb = structts2vec(shapes, {'zb'}, t(i));
  [~,p] = jphi_estimate(rb, zb, ip, tok.rg, tok.zg, opts);
  pcurrt.Data(i,:) = p(:);
end