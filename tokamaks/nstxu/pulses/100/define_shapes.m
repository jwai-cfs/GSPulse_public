function shapes = define_shapes(settings, tok)
%
%  inputs: opts struct
%          if opts.plotlevel >= 1, makes a gui plot of the target shapes
%          if opts.plotlevel >= 2, also plots the shape target timetraces
% 
%  outputs: the shape struct that has the target shapes vs time 
%
%     (rb.Data, rb.Time) - R position of boundary target shape, 
%                          timebase for R position of boundary target shape
%
%     (zb.Data, zb.Time)         - Z position of boundary target shape
%     (rx.Data, rx.Time)         - R position of target x-point
%     (zx.Data, zx.Time)         - Z position of target x-point
%     (rtouch.Data, rtouch.Time) - R position of target touch point
%     (ztouch.Data, ztouch.Time) - Z position of target touch point
%     (rbdef.Data, rbdef.Time)   - R of point where boundary flux (psibry)
%                                    is evaluated
%     (zbdef.Data, zbdef.Time)   - Z of " " "
%
%     Note that each of these can use different timebases which are later
%     interpolated. Each quantity should be defined at each time, but will 
%     only enter the optimization depending on the optimization weights. 


% some reference equilibria to draw new shapes from
eqdir = [getenv('GSROOT') '/tokamaks/nstxu/data/'];
eqLIM = load([eqdir 'eqLIM']).eq;  
eqLSN = load([eqdir 'eqLSN']).eq;


% Define target shapes
i = 0;
rb = {};
zb = {};
rmin = min(tok.limdata(2,:));


% start with a small limited plasma at t=0
eq = eqLIM;
i = i + 1;
t(i) = 0;
s = shape_analysis(eq.rbbbs, eq.zbbbs);
s.zsurf = -0.02;
s.aminor = 0.48;
s.triu = 0.2;
s.tril = 0.2;
[rb{i}, zb{i}] = shape_edit(eq.rbbbs, eq.zbbbs, s);
rb{i} = rb{i} - min(rb{i}) + rmin;


% keep same small limited plasma at t=0.1
i = i + 1;
t(i) = 0.1;
rb{i} = rb{i-1};
zb{i} = zb{i-1};


% plasma is about to divert at t=0.2
% larger limited plasma but with an x-point
eq = eqLSN;
i = i + 1;
t(i) = 0.2;
s = shape_analysis(eq.rbbbs, eq.zbbbs);
s.zsurf = -0.02;
s.elong = 1.6;
s.triu = 0.3;
s.squo = -0.05;
s.squi = -0.0;
s.sqli = -0.31;
s.sqlo = -0.2;
[rb{i}, zb{i}] = shape_edit(eq.rbbbs, eq.zbbbs, s);
rb{i} = rb{i} - min(rb{i}) + rmin;

% divert plasma by moving radially from wall at t=0.23
i = i + 1;
t(i) = 0.29;
rb{i} = rb{i-1} + 0.03;
zb{i} = zb{i-1};


% keep the same shape until shot end at t=1
i = i + 1;
t(i) = 1;
rb{i} = rb{i-1};
zb{i} = zb{i-1};
t = t(:);

%% Map target shapes to boundary-control points
% At this point, each of the (rb{i}, zb{i}) define a shape. Now sort these
% and map them onto control points. 


% interpolate to finer boundary
for i = 1:length(rb)
  [rb{i}, zb{i}] = interparc(rb{i}, zb{i}, 200, 1, 0);  % interpolate
  [rb{i}, zb{i}] = sort_ccw(rb{i}, zb{i}, 0.9, 0);      % sort 
end
  

% define control segments
segopts.rc = 0.85; 
segopts.zc = 0;
segopts.a = 0.3;
segopts.b = 0.5;
segopts.plotit = 0;
segopts.seglength = 4;
segs = gensegs(40, segopts);


% find intersections of boundary with segments
rcp = [];
zcp = [];
ncp = length(rb);
for i = 1:ncp
  [rcp(i,:), zcp(i,:)] = seg_intersections(rb{i}, zb{i}, segs, 0);  
end

shapes.rb.Time = t;
shapes.rb.Data = rcp;

shapes.zb.Time = t;
shapes.zb.Data = zcp;


%% define x-points, touch points
% for time>0.2, just take the lowest boundary z-position as the x-point
rx = [];
zx = [];

for i = 1:length(rb)
  if t(i) >= 0.2
    [zx(i), j] = min(zb{i});
    rx(i) = rb{i}(j);
  else
    rx(i) = nan;
    zx(i) = nan;
  end
end
rx = inpaint_nans(rx(:));
zx = inpaint_nans(zx(:));

shapes.rx.Time = t;
shapes.rx.Data = rx;

shapes.zx.Time = t;
shapes.zx.Data = zx;

shapes.rtouch.Time = [0 1]';
shapes.rtouch.Data = [0.315 0.315]';

shapes.ztouch.Time = [0 1]';
shapes.ztouch.Data = [0 0]';

shapes.rbdef = shapes.rb;
shapes.zbdef = shapes.zb;

shapes.rstrike.Time = [0 1]';
shapes.rstrike.Data = [nan nan]';

shapes.zstrike.Time = [0 1]';
shapes.zstrike.Data = [nan nan]';

% shapes.rbdef.Time = t;
% shapes.rbdef.Data = rx(:);
% shapes.rbdef.Data(t<=0.2) = 0.315;
% 
% shapes.zbdef.Time = t;
% shapes.zbdef.Data = zx(:);
% shapes.zbdef.Data(t<=0.2) = 0;

% plots
if settings.plotlevel >= 2
  summary_shape_plot(shapes, tok);
  plot_structts(shapes, fieldnames(shapes), [], 3);
  sgtitle('Shape targets'); drawnow
end


















