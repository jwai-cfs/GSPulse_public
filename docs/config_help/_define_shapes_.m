% USAGE:    
%   shapes = define_shapes(tok)
%
% DESCRIPTION: 
%
% Define the time-dependent shaping parameters. 
%
% INPUTS: 
%  tok - the tokamak geometry object
%
% OUTPUTS: 
%  shapes - struct with fields defining the time-evolution of the shapes.
%           See below for details. 
%
% ADDITIONAL INFO
%   
% GSPulse contains a GUI for defining and editing shapes. To run the GUI, 
% in a terminal: 
%
%    cd to the GSPulse root directory
%    python shape_gui/shape_gui.py -tok <tokamak> 
%
% At this point, the supported tokamaks for the GUI are SPARC, ARC, MAST-U,
% and NSTX-U. One can use the GUI to define the shapes and save them to
% file. 
% 
% In this script we define the time evolution of the shapes. To do this,
% define a timebase, and a sequence of shapes corresponding to this
% timebase. Then, do
%
% shapes = shapefiles2shapes(shapefiles, tok, time)
%
% This workflow will populate the relevant fields of the shapes object. The
% fields in the shapes struct are: 
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
%     Each of these can use different timebases which are later
%     interpolated. Each quantity should be defined at each time, but will 
%     only enter the optimization depending on the optimization weights. 


