function tok = define_tok(~)

% load geometry
fn = [getenv('GSROOT') '/tokamaks/nstxu/device/nstxu_obj_config.mat'];
tok_data_struct = load(fn).tok_data_struct;

fn = [getenv('GSROOT') '/tokamaks/nstxu/device/sysid_fits.mat'];
sysid_fits = load(fn).sysid_fits;

% connect tok_data_struct into circuits instead of individual coils
addpath([getenv('GSROOT') '/tokamaks/nstxu/device'])
circ = nstxu2016_circ(tok_data_struct);

% these are mutual inductances and resistances that were obtained from
% fitting data from several shots. Overwrite the toksys-calculated values
% with these fitted values. 
mxx = sysid_fits.Mxx;
rxx = sysid_fits.Rxx;

iremove = [3 4 7 11 12];
mxx(iremove,:) = [];
mxx(:,iremove) = [];
rxx(iremove) = [];

tok_data_struct.mvv     = mxx(circ.iivx, circ.iivx);
tok_data_struct.mcc     = mxx(circ.iicx, circ.iicx);
tok_data_struct.mcv     = mxx(circ.iicx, circ.iivx);
tok_data_struct.resv    = rxx(circ.iivx);
tok_data_struct.resc    = rxx(circ.iicx);
tok_data_struct.mpc     = tok_data_struct.mpc * circ.Pcc;
tok_data_struct.mpv     = tok_data_struct.mpv * circ.Pvv;
tok_data_struct.ccnames = circ.ccnames;
tok_data_struct.nc      = length(circ.iicx);
tok_data_struct.nv      = length(circ.iivx);

% map from "tok_data_struct" object (toksys format) to "tok" object
% (gspulse format) which is mostly the same but with a few stripped fields
tok = tok_data_struct2tok(tok_data_struct);
