function segs = gensegs(n, opts)
% =========================================================================
% Description: 
% generate equidistant shape control segments
%
% Method:
% makes a small ellipse defined by (opts.a, opts.b) centered at 
% (opts.rc, opts.zc). Makes a larger ellipse that is opts.seglength
% bigger. Connects these two ellipses to form the segments
% 
% Example:
%   opts.plotit = 1;
%   segs = gensegs(40, opts)
% =========================================================================

if ~exist('opts','var'), opts = struct; end
if ~isfield(opts, 'rc'), opts.rc = 0.85; end
if ~isfield(opts, 'zc'), opts.zc = 0; end
if ~isfield(opts, 'a'),  opts.a = 0.3; end
if ~isfield(opts, 'b'),  opts.b = 0.5; end
if ~isfield(opts, 'seglength'), opts.seglength = 4; end
if ~isfield(opts, 'plotit'), opts.plotit = false; end

th = linspace(0, 2*pi, 200)';

rin = opts.rc + opts.a*cos(th);
zin = opts.zc + opts.b*sin(th);

rout = opts.rc + opts.seglength * opts.a * cos(th);
zout = opts.zc + opts.seglength * opts.b * sin(th);

[r0, z0] = interparc(rout, zout, n, 0, 0);

k = dsearchn([rin zin], [r0 z0]);
rf = rin(k);
zf = zin(k);

segs = [r0 rf z0 zf];

if opts.plotit
  plot([r0 rf]', [z0 zf]', 'k')
  axis equal  
  hold on
  i = 1;
  plot([r0(i) rf(i)]', [z0(i) zf(i)]', 'g', 'linewidth', 2)
end