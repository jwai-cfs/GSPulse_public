% Create the connected tok_data_struct object from the external data
% sources for nstxu. 
clear; clc; close all

saveit = 0;
savefn = [getenv('GSROOT') '/tokamaks/nstxu/tok/nstxu_tok.mat'];

% load geometry
vac_sys = load('NSTXU_vacuum_system.mat').NSTXU_vacuum_system;
tok_data_struct = vac_sys.build_inputs.tok_data_struct;
tok_data_struct.imks = 1;
circ = nstxu2016_circ(tok_data_struct);


% mutual inductances and resistances
mxx = vac_sys.sysid_fits.Mxx;
rxx = vac_sys.sysid_fits.Rxx;
mvv = mxx(circ.iivx, circ.iivx);
mcc = mxx(circ.iicx, circ.iicx);
mcv = mxx(circ.iicx, circ.iivx);
resv = rxx(circ.iivx);
resc = rxx(circ.iicx);
mpc = tok_data_struct.mpc * circ.Pcc;
mpv = tok_data_struct.mpv * circ.Pvv;
% mpp = (tok_data_struct.mpp + tok_data_struct.mpp') / 2;
mpp = load('nstxu_obj_config2016_6565.mat').tok_data_struct.mpp;


ccnames = circ.ccnames(:);
nc = circ.ncx;
nv = circ.nvx;
rg = tok_data_struct.rg;
zg = tok_data_struct.zg;
nr = tok_data_struct.nr;
nz = tok_data_struct.nz;
limdata = tok_data_struct.limdata;

% descriptions
d.mcc = 'mutual inductances from coils to coils [Wb/A]';
d.mcv = 'mutual inductances from coils to vessels [Wb/A]';
d.mvv = 'mutual inductances from vessels to vessels [Wb/A]';
d.mpc = 'mutual inductances from plasma grid to coils [Wb/A]';
d.mpv = 'mutual inductances from plasma grid to vessels [Wb/A]';
d.mpp = 'mutual inductances from plasma grid to plasma grid [Wb/A]';
d.resc = 'coil resistances [Ohms]';
d.resv = 'vessel resistances [Ohms]';
d.ccnames = 'coil names';
d.nc = 'number of coils';
d.nv = 'number of vessel elements';
d.rg = 'radial positions on grid [m]';
d.zg = 'vertical positions on grid [m]';
d.nr = 'grid size (radial)';
d.nz = 'grid size (vertical)';
d.limdata = '(z,r) of vertices defining limiter';
descriptions = d;

% save data
tok = variables2struct(mcc, mcv, mvv, mpc, mpv, mpp, resc, resv, ...
  ccnames, nc, nv, rg, zg, nr, nz, limdata, descriptions); 

if saveit
  save(savefn, 'tok')
end
































