function [pprime, ttprime, press] = normalize_profile_shape(...
  npsin, pprime, ttprime, press_in)
% *************************************************************************
% Description: 
% Normalize profile shapes to a consistent format. In general, FBT doesn't
% care if the profile shapes are not normalized (e.g. scaled) since it just
% scales the corresponding coefficients. 
%
% However, we do need to normalize in some cases because the scaling can
% cause issues. For example, if using different sources of profile shapes
% (GSPulse and RAPTOR) and taking an average, this is done incorrectly if
% the profiles are not normalized first. 
% *************************************************************************

% ensure all profiles have length npsin
pprime = interp1(linspace(0,1,length(pprime)), pprime, linspace(0,1,npsin)); 
ttprime = interp1(linspace(0,1,length(ttprime)), ttprime, linspace(0,1,npsin)); 

% pressure basis from integrating pprime 
press_from_pprime = bfprmex(pprime(:) * ones(1,3)) * [1 0 0]';  

if ~exist('press_in', 'var') || isempty(press_in)  
  press = press_from_pprime;
else
  % some scale and shift transformations for setting finite edge pressure
  press_in = interp1(linspace(0,1,length(press_in)), press_in, linspace(0,1,npsin)); 
  press_in_zero_edge = press_in - press_in(end);
  scale = press_in_zero_edge(:) \ press_from_pprime(:);
  press = press_in * scale;
end

% normalize to have unit length and positive
scale = 1 / norm(pprime) * median(sign(pprime));
pprime = pprime * scale;   
press = press * scale;
ttprime = ttprime / norm(ttprime) * median(sign(ttprime));

  