function shapes = shapefiles2shapes(shapefns, tok, t)
% =========================================================================
% Description: 
%   import the shapes created by the python shape GUI onto waveforms for
%   the specified time basis
%
% Inputs: 
%  shapefns - cellarray with the shape filenames to load
%  tok           - tokamak geometry struct, see help _define_tok.m
%  t  - time basis corresponding to each shape
%
% Outputs: 
%  shapes - struct of timeseries, containing a timeseries for all of the
%    various shaping parameters
%
% =========================================================================
t = t(:);

% load data
for i = 1:length(shapefns)
  try
    [~, ~, rcp, zcp, rstrike, zstrike, rx, zx] = load_shape_file(shapefns{i}, 0);
  catch
    warning('Could not load %s', shapefns{i});
  end

  shapes.rb.Data(i,:) = rcp;
  shapes.zb.Data(i,:) = zcp;
  shapes.rstrike.Data(i,:) = rstrike;
  shapes.zstrike.Data(i,:) = zstrike;
  shapes.rx.Data(i,:) = rx;
  shapes.zx.Data(i,:) = zx;   
end
shapes.rbdef = shapes.rb;
shapes.zbdef = shapes.zb;


% assign timebase to shapes
shapes.rb.Time = t;
shapes.zb.Time = t;
shapes.rstrike.Time = t;
shapes.zstrike.Time = t;
shapes.rx.Time = t;
shapes.zx.Time = t;
shapes.rtouch.Time = t;
shapes.ztouch.Time = t;
shapes.rbdef.Time = t;
shapes.zbdef.Time = t;

% assign touch point at z=0 on inboard wall limiter
riwl = seg_intersections(tok.rl, tok.zl, [min(tok.rl)-0.1 min(tok.rl)+sqrt(eps) 0 0]);
shapes.rtouch.Data = ones(size(t)) * riwl;
shapes.ztouch.Data = ones(size(t)) * 0;

% transpose data if necessary
shapes = check_structts_dims(shapes);