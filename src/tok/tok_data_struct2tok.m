function tok = tok_data_struct2tok(tds)
% =========================================================================
% Description: 
%  convert from the TokSys-formatted "tok_data_struct" object to the "tok"
%  object, which is what GSPulse uses internally for geometry and mutual
%  inductance definitions. 
%
% 
% Inputs: 
%  tds - tok_data_struct object from the TokSys code
%
% Outputs:
%  tok - tokamak geometry object used by GSPulse. 
%        See also, help _define_tok.m
%
% =========================================================================

% copy a bunch of fields directly over with no changes
fds2copy = {'rg', 'zg', 'nr', 'nz', 'nc', 'nv', 'limdata', 'mcc', ...
  'mcv', 'mvv', 'mpc', 'mpv', 'resc', 'resv'};
tok = copyfields(struct, tds, fds2copy, 0);

% copy some fields with modifications
tok.ccnames = cellstr(tds.ccnames);
tok.rl = tds.limdata(2,:)';
tok.zl = tds.limdata(1,:)';

% mpp is stored either in full format [nz*nr x nz*nr] or compressed format
% [nz*nr x nr]. Unwrap it to full format if compressed.
if size(tds.mpp,2) ~= tds.nr*tds.nz
  tok.mpp = unwrap_mpp(tds.mpp, tds.nz, tds.nr);
else
  tok.mpp = tds.mpp;
end

[tok.rgg, tok.zgg] = meshgrid(tok.rg, tok.zg);
in = inpolygon(tok.rgg, tok.zgg, tok.rl, tok.zl);
tok.outside_vessel_mask = ~in;