function plot_coil_labels(tok, names, cccirc, dr, dz)
% =========================================================================
% Description: 
%  plot coil labels on a figure
%
% Inputs: 
%  tok - tokamak geometry object
%  names - cell array of coil names
%  ccccirc - circuit connection vector
%  (dr,dz) - how much to shift label text
%  
% Outputs: 
%   None
%
% =========================================================================

if ~exist('dr', 'var'), dr = 0; end
if ~exist('dz', 'var'), dz = 0; end
if ~exist('cccirc', 'var'), cccirc = [1 1 2 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 -19]; end

fcdata = tok.fcdata;
ecdata = [];
if isfield(tok, 'ecdata') && ~isempty(tok.ecdata)
  ecdata = tok.ecdata;
  xmin = min(ecdata(2,:) - ecdata(4,:)/2);
  xmax = max(ecdata(2,:) + ecdata(4,:)/2);
  ymin = min(ecdata(1,:) - ecdata(3,:)/2);
  ymax = max(ecdata(1,:) + ecdata(3,:)/2);
  xcenter = (xmin + xmax)/2;
  ycenter = (ymin + ymax)/2;
  dx = xmax - xmin;
  dy = ymax - ymin;
  ecdata = [ycenter xcenter dy dx 0 0]';  
end
ccdata = [ecdata fcdata];
nc = size(ccdata,2);

for i = 1:length(names)
  
  k = find(cccirc == i | cccirc == -i);
  z = ccdata(1,k);
  r = ccdata(2,k);

  zcu = mean(z(z>=0));
  zcl = mean(z(z<0));
  rc = mean(r);
  
  rc = rc + dr;
  zcu = zcu + dz;
  zcl = zcl + dz;

  t = text(rc, zcu, names{i});
  t.FontSize = 12;
  t.FontWeight = 'bold';
  t.HorizontalAlignment = 'center';
  t.BackgroundColor = [1 1 1 0.3];


  t = text(rc, zcl, names{i});
  t.FontSize = 12;
  t.FontWeight = 'bold';
  t.HorizontalAlignment = 'center';
  t.BackgroundColor = [1 1 1 0.3];

end