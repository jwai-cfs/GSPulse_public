function h = plot_efit_region(z, r, dz, dr, ac, ac2, rgb, varargin)
%
% PLOT_EFIT_REGION Function for plotting a parallelogram with geometry
%                  described by the EFIT conductor representation.
%
% USAGE: plot_efit_region.m
%
% INPUTS:
%
%   z..........vertical position of parallelogram center [m]
%   r..........major radius of parallelogram center      [m]
%   dz.........full height of the parallelogram          [m]    
%   dr.........full width of the parallelogram           [m]
%   ac.........counterclockwise rotation (angled bottom) [deg]
%   ac2........counterclockwise rotation (flat bottom)   [deg]
%   rgb........rgb color vector
%   varargin ... additional arguments to pass to plotting function (fill)
%
% OUTPUTS: 
%
%   NONE
%
%..........................................................................

% Convert angles to radians
ac  = ac  * pi/180;
ac2 = ac2 * pi/180;

% The parallelogram is a rectangle
if ac == 0 && ac2 == 0
    r1 = r - dr/2;
    r2 = r + dr/2;
    r3 = r + dr/2;
    r4 = r - dr/2;
    
    z1 = z - dz/2;
    z2 = z - dz/2;
    z3 = z + dz/2;
    z4 = z + dz/2;   
end

% The parallelogram has an angled bottom
if ac ~= 0 && ac2 == 0 
    r1 = r - dr/2;
    r2 = r + dr/2;
    r3 = r + dr/2;
    r4 = r - dr/2;
        
    z1 = z - sign(ac)*dz/2 - dr*tan(ac)/2;
    z2 = z - sign(ac)*dz/2 + dr*tan(ac)/2;
    z3 = z + sign(ac)*dz/2 + dr*tan(ac)/2;
    z4 = z + sign(ac)*dz/2 - dr*tan(ac)/2; 
end
    
% The parallelogram has a flat bottom
if ac == 0 && ac2 ~= 0
    r1 = r - dz/tan(ac2)/2 - dr/2;
    r2 = r - dz/tan(ac2)/2 + dr/2;
    r3 = r + dz/tan(ac2)/2 + dr/2;
    r4 = r + dz/tan(ac2)/2 - dr/2;

    z1 = z - dz/2;
    z2 = z - dz/2;
    z3 = z + dz/2;
    z4 = z + dz/2;
    
end

% Plot the parallelogram
h = fill([r1 r2 r3 r4],  [z1 z2 z3 z4], rgb, varargin{:});

end
