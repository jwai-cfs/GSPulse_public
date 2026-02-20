function [ax,ay,az]= dir_cosine(x1,y1,z1,x2,y2,z2)
% =========================================================================
% Description: 
%  find the direction cosines (contributions of a component of the basis 
%  to a unit vector in that direction) for the vector that starts at 
%  (x1,y1,z1) and ends at (x2,y2,z2)
%
% Inputs: 
% (x1,y1,z1) - starting coordinate
% (x2,y2,z2) - ending coordinate
% 
% Outputs: 
%  (ax,ay,ax) - direction cosines for the 3 principal directions
%  
% =========================================================================
d= sqrt((x2-x1).^2+(y2-y1).^2+(z2-z1).^2);
ax= (x2-x1)./d;
ay= (y2-y1)./d;
az= (z2-z1)./d;