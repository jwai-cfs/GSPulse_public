function cmats = output_model(dpsizrdx, tok, shapes, t, fds2control)
% =========================================================================
% Description:
%
%
%
% Inputs:
%
%
%
% Outputs:
%
%
%
% Additional info:
%
%
% =========================================================================
N = length(t);
cmats = cell(N,1);

for i = 1:N
  ref = structts2struct(shapes, fieldnames(shapes), t(i));
  cdata = build_cmat(dpsizrdx, ref, tok, fds2control);
  cmats{i} = struct2vec(cdata, fds2control);
end
end


function c = build_cmat(dpsizrdx, ref, tok, fds2control)

c = struct;   % this will hold all the derivatives

% response of currents is identity
if ismember('ic', fds2control)
  c.ic = [eye(tok.nc) zeros(tok.nc, tok.nv)];
end
if ismember('iv', fds2control)
  c.iv = [zeros(tok.nv,tok.nc) eye(tok.nv)];
end

% linear combinations of currents
if ismember('Aic', fds2control)
  c.Aic = [tok.Aic zeros(size(tok.Aic,1), tok.nv)];
end

% isoflux response at control points
if any(contains(fds2control, 'psicp'))
  c.psicp = multigrid2pt(tok.rg, tok.zg, dpsizrdx, ref.rb, ref.zb);
end

% isoflux response at strike points (note: using the term "strike point"
% loosely, could be any point not touching the plasma)
if any(contains(fds2control, 'psisp'))
  c.psisp = multigrid2pt(tok.rg, tok.zg, dpsizrdx, ref.rstrike, ref.zstrike);
end

% x-point responses
if ismember('psix_r', fds2control) || ismember('psix_z', fds2control)
  [c.psix, c.psix_r, c.psix_z] = multigrid2pt(tok.rg, tok.zg, dpsizrdx, ref.rx, ref.zx);
end

% touch point response
if any(contains(fds2control, 'psitouch'))
  c.psitouch = multigrid2pt(tok.rg, tok.zg, dpsizrdx, ref.rtouch, ref.ztouch);
end

% touch point response
if any(contains(fds2control, 'psibry'))
  c.psibry = mean(multigrid2pt(tok.rg, tok.zg, dpsizrdx, ref.rbdef, ref.zbdef));
end

% isoflux differences
if ismember('diff_psicp_psix1', fds2control)
  c.diff_psicp_psix1 = c.psicp - c.psix(1,:);
end
if ismember('diff_psicp_psix2', fds2control)
  c.diff_psicp_psix2 = c.psicp - c.psix(2,:);
end
if ismember('diff_psicp_psitouch', fds2control)
  c.diff_psicp_psitouch = c.psicp - c.psitouch;
end
if ismember('diff_psicp_psibry', fds2control)
  c.diff_psicp_psibry = c.psicp - c.psibry;
end
if ismember('diff_psisp_psix1', fds2control)
  c.diff_psisp_psix1 = c.psisp - c.psix(1,:);
end
if ismember('diff_psisp_psix2', fds2control)
  c.diff_psisp_psix2 = c.psisp - c.psix(2,:);
end
if ismember('diff_psix1_psix2', fds2control)
  c.diff_psix1_psix2 = c.psix(1,:) - c.psix(2,:);
end
if ismember('diff_psix2_psix3', fds2control)
  c.diff_psix2_psix3 = c.psix(2,:) - c.psix(3,:);
end
end