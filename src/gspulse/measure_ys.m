function y = measure_ys(psizr, ic, fds2control, shapes, tok, t)
% =========================================================================
% Description:
%  measure all of the controlled outputs (specified by fds2control) for
%  a given flux distribution
%
% Inputs:
%  psizr - array with size [nz*nr x length(t)] describing the flux on the
%          grid at each time
%  ic    - coil currents with size [ncoils x length(t)] describing the coil
%          currents vs time
%  fds2control - same as settings.fds2control, see help _define_settings.m
%  shapes        - shapes struct, see help _define_shapes.m
%  tok           - tokamak geometry struct, see help _define_tok.m
%  t - timebase for the inputs and outputs
%
% Outputs:
%  y  - struct of timeseries, with waveforms defined for each of the
%  fds2control
%
% =========================================================================
N = length(t);
y = cell(N,1);
for i = 1:N
  ref = structts2struct(shapes, fieldnames(shapes), t(i));
  ydata = measure_y_fun(psizr(:,i), ref, tok, ic(:,i), fds2control);
  ydata = check_struct_dims(ydata);
  y{i} = struct2vec(ydata, fds2control);
end

end

function y = measure_y_fun(psizr, ref, tok, ic, fds2control)

y = struct;   % this will hold all the signals

% measurement of currents is identity
if ismember('ic', fds2control)
  y.ic = ic;
end

% linear combinations of currents
if ismember('Aic', fds2control)
  y.Aic = tok.Aic * ic;
end

% isoflux measurement at control points
if any(contains(fds2control, 'psicp'))
  y.psicp = bicubicHermite(tok.rg, tok.zg, psizr, ref.rb, ref.zb);
end

% isoflux measurement at strike points (note: using the term "strike point"
% loosely, could be any point not touching the plasma)
if any(contains(fds2control, 'psisp'))
  y.psisp = bicubicHermite(tok.rg, tok.zg, psizr, ref.rstrike, ref.zstrike);
end

% x-point responses
if ismember('psix', fds2control) || ismember('psix_r', fds2control) || ismember('psix_z', fds2control)
  [y.psix, y.psix_r, y.psix_z] = bicubicHermite(tok.rg, tok.zg, psizr, ref.rx, ref.zx);
end

% touch point measurement
if any(contains(fds2control, 'psitouch'))
  y.psitouch = bicubicHermite(tok.rg, tok.zg, psizr, ref.rtouch, ref.ztouch);
end

% touch point measurement
if any(contains(fds2control, 'psibry'))
  y.psibry = mean(bicubicHermite(tok.rg, tok.zg, psizr, ref.rbdef, ref.zbdef));
end

% isoflux differences
if ismember('diff_psicp_psix1', fds2control)
  y.diff_psicp_psix1 = y.psicp - y.psix(1);
end
if ismember('diff_psicp_psix2', fds2control)
  y.diff_psicp_psix2 = y.psicp - y.psix(2);
end
if ismember('diff_psicp_psitouch', fds2control)
  y.diff_psicp_psitouch = y.psicp - y.psitouch;
end
if ismember('diff_psicp_psibry', fds2control)
  y.diff_psicp_psibry = y.psicp - y.psibry;
end
if ismember('diff_psisp_psix1', fds2control)
  y.diff_psisp_psix1 = y.psisp - y.psix(1);
end
if ismember('diff_psisp_psix2', fds2control)
  y.diff_psisp_psix2 = y.psisp - y.psix(2);
end
if ismember('diff_psix1_psix2', fds2control)
  y.diff_psix1_psix2 = y.psix(1) - y.psix(2);
end
if ismember('diff_psix2_psix3', fds2control)
  y.diff_psix2_psix3 = y.psix(2) - y.psix(3);
end
end