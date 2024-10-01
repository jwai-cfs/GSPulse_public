
% convert from MEQ naming schema of geometry (L struct) to GSPulse naming
% schema of geometry (tok struct)

function tok = L2tok(L)

% convert from L to tok
tok = struct;
tok.ccnames = L.G.dima;
tok.limdata = [L.G.zl L.G.rl]';
tok.rl = L.G.rl;
tok.zl = L.G.zl;
tok.mcc = L.G.Maa;
tok.mcv = L.G.Mau;
tok.mpc = L.G.Mxa;
tok.mpp = L.G.Mxx;
tok.mpv = L.G.Mxu;
tok.mvv = L.G.Muu;
tok.nr = length(L.G.rx);
tok.nz = length(L.G.zx);
tok.nc = L.G.na;
tok.nv = L.G.nu;
tok.resc = L.G.Ra;
tok.resv = L.G.Ru;
tok.rg = L.G.rx;
tok.zg = L.G.zx;
tok.descriptions = '';
[tok.rgg, tok.zgg] = meshgrid(tok.rg, tok.zg);
in = inpolygon(tok.rgg, tok.zgg, tok.limdata(2,:), tok.limdata(1,:));
tok.outside_vessel_mask = ~in;