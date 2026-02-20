function Pcc = cccirc_to_Pcc(cccirc)
% =========================================================================
% Description: 
%  make a circuit connection matrix from cccirc. The output matrix is such
%  that:
%
%  ic = Pcc * icx
%  icx = pinv(Pcc) * ic
%  
%  where icx is the circuit currents and ic is the individual coil currents
% 
% Inputs: 
%  cccirc - circuit connection vector. As an example,
% 
%     cccirc = [1 1 1 2 -2 3]
% 
%  indicates that the first 3 coils are wired in series, the next 2 coils
%  are wired in anti-series, and the last coil is independent. 
%
% =========================================================================
n = max(abs(cccirc));
Pcc = zeros(length(cccirc), n);  

for i = 1:n
  Pcc(cccirc == i,i) = 1;
  Pcc(cccirc == -i,i) = -1;  
end