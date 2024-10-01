function plot_nstxu_geo(tok_data_struct, options)
%
% PLOT_NSTXU_GEO
%
% SYNTAX: plot_nstxu_geo.m
% 
% PURPOSE: Plot NSTX-U geometry.
% 
% INPUTS:
%
%   tok_data_struct....TokSys struct describing the tokamak geometry
%
%   options............structure containing any of the following optional
%                      fields:
%
%       iaxes.....show axes and labels: no [0] or yes [1]
%
%       icolor....color conductor cross-sections based upon material:
%                 no color [0] or color [1]
%
%       igrid.....plot default grid: no [0] or yes [1]
%
%       ifl.......plot the flux loops: no [0] or yes [1]
%
%       ifull.....plot the full cross-section (left + right halves): 
%                 no [0] or yes [1]
%
% OUTPUTS: Plot of NSTX-U geometry with formatting specified by options.
%..........................................................................
 
% Unpack tok_data_struct and options

ecdata  = tok_data_struct.ecdata;
fcdata  = tok_data_struct.fcdata;
vvdata  = tok_data_struct.vvdata;
limdata = tok_data_struct.limdata;
fcnames = tok_data_struct.fcnames;

% etav = tok_data_struct.make_tok_inputs.etav;
if nargin == 1, options = struct; end
if ~isfield(options, 'iaxes'),  options.iaxes  = 1; end
if ~isfield(options, 'icolor'), options.icolor = 1; end
if ~isfield(options, 'igrid'),  options.igrid  = 0; end
if ~isfield(options, 'ifl'),    options.ifl    = 0; end
if ~isfield(options, 'ifull'),  options.ifull  = 0; end
if ~isfield(options, 'plot_vv'),  options.plot_vv = 1; end
if ~isfield(options, 'labels'),  options.labels = 0; end
 
% scrsz = get(groot,'ScreenSize');
% figure('Position',[1 scrsz(4)/2 scrsz(3)/4 scrsz(4)])

%.................
% Plot the OH coil

% Just plot the shape of the OH coil and not each individual turn

zbotleft  = min(ecdata(1,:)) - ecdata(3,1)/2;
rbotleft  = min(ecdata(2,:)) - ecdata(4,1)/2;
ztopleft  = max(ecdata(1,:)) + ecdata(3,1)/2;
rbotright = max(ecdata(2,:)) + ecdata(4,1)/2;

dz = ztopleft  - zbotleft;
dr = rbotright - rbotleft;

position = [rbotleft zbotleft dr dz];

if options.icolor
    rgb = [1 .75 0]; % Copper(OFHC)
else
    rgb = [1 1 1]*0.9;
end

rectangle('Position', position, 'FaceColor', rgb)

if options.ifull
    
    zbotleft  =  min(ecdata(1,:)) - ecdata(3,1)/2;
    rbotleft  = -max(ecdata(2,:)) - ecdata(4,1)/2;
    ztopleft  =  max(ecdata(1,:)) + ecdata(3,1)/2;
    rbotright = -min(ecdata(2,:)) + ecdata(4,1)/2;
    
    dz = ztopleft  - zbotleft;
    dr = rbotright - rbotleft;
    
    position = [rbotleft zbotleft dr dz];
    rectangle('Position', position, 'FaceColor', rgb)
    
end

%..................
% Plot the PF coils

for ii = 1:size(fcdata,2)
    
    zbotleft = fcdata(1,ii) - fcdata(3,ii)/2;
    rbotleft = fcdata(2,ii) - fcdata(4,ii)/2;
    position = [rbotleft zbotleft fcdata(4,ii) fcdata(3,ii)];
    rectangle('Position', position, 'FaceColor', rgb) 
    
    if options.labels
      text(rbotleft, fcdata(1,ii), fcnames(ii,:), 'fontsize', 10, 'fontweight', 'bold');
    end
    
    if options.ifull
        
        zbotleft =  fcdata(1,ii) - fcdata(3,ii)/2;
        rbotleft = -fcdata(2,ii) - fcdata(4,ii)/2;
        position = [rbotleft zbotleft fcdata(4,ii) fcdata(3,ii)];
        rectangle('Position', position, 'FaceColor', rgb)
        
    end
        
end
    
%............................
% Plot the passive conductors 

z   = vvdata(1,:);
r   = vvdata(2,:);
dz  = vvdata(3,:);
dr  = vvdata(4,:);
ac  = vvdata(5,:);
ac2 = vvdata(6,:);

hold on
for ii = 1:size(vvdata,2)
    
%     if options.icolor
%         switch etav(ii)
%             case 0.0209
%                 rgb = [1 0.75 0]; 
%             case 1.2800
%                 rgb = [1 0 0];
%             case 0.7200
%                 rgb = [0 0 1];
%             case 0.7400
%                 rgb = [0 0 1];
%         end
%     else
%         rgb = [1 1 1];
%     end
    if options.plot_vv
      plot_efit_region(z(ii), r(ii), dz(ii), dr(ii), ac(ii), ac2(ii), rgb);
      hold on

      if options.ifull
          plot_efit_region(z(ii), -r(ii), dz(ii), dr(ii), -ac(ii), ...
              -ac2(ii), rgb);
          hold on
      end
    end

end

%.................
% Plot the limiter

zlim = limdata(1,:);
rlim = limdata(2,:);


plot(rlim, zlim, 'color', [1 1 1]*0.4, 'LineWidth', 3)
hold on

if options.ifull
    hold on
    plot(-rlim, zlim, '-k', 'LineWidth', 2)
end

%................
% Plot flux loops

if options.ifl
    z = tok_data_struct.fldata(1,:);
    r = tok_data_struct.fldata(2,:);
    
    hold on
    scatter(r, z, 'or', 'filled')
end

%................
% Figure settings

axis equal

if options.ifull
    axis([-2.2 2.2 -2.2 2.2])
else
    axis([0.1 2.1 -2 2])
end 

if options.iaxes
    hXl = xlabel('r [m]');
    hYl = ylabel('z [m]');
    hTl = title('NSTX-U');
    set([hXl, hYl], 'FontSize', 14, 'FontWeight', 'bold')
    set(hTl, 'FontSize', 14.5, 'FontWeight', 'bold') 
else
    axis off
end

if options.igrid
    grid on
end

end
      