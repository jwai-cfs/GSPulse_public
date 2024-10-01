% Form the grouping circuits for the coils and vacuum vessel elements. 
% 
% cccirc is coil circuit connections. For example: [1 -1 2 2 2 3] says that
% coils 1 and 2 are connected in antiseries, coils 3-5 are connected in
% series, and coil 6 is independent. 
%
% vvgroup and vvcirc are similar, but subtly different. vvgroup = [1 1 2 3]
% says that elements 1 and 2 are lumped together (in parallel) into
% group 1, elements 3 = group 2, and element 4 = group 3. vvcirc = [1 2 2]
% now applies on the groups, and says that group 1 is independent and
% groups 2 and 3 are connected in series. However, there are almost never 
% series connections in the vessel groupings ==> vvcirc = 1:max(vvgroup); 
%
% Pcc, Pvv, and Pxx are transition maps from connected to unconnected
% circuits. Let ic, iv be unconnected coil and vessel current vectors. 
% Let icx, ivx be connected coil circuit and vessel circuit vectors. Let 
% is := [ic; iv; ip] with ip=plasma current, isx = [icx; ivx; ip]. 
%
% Then:  ic  = Pcc * icx
%        icx = pinv(Pcc) * ic
%        iv  = Pvv * ivx
%        ivx = pinv(Pvv) * iv
%
% Josiah Wai

function circ = nstxu2016_circ(tok_data_struct)


cccirc = [1 2 0 0 3 3 4 4 4 4 0 0 0 0 0 0 5 5 5 5 6 6 6 6 7 7 0 0 8];
Pcc = cccirc_to_Pcc(cccirc);
ccnames = {'OH', 'PF1AU', 'PF2U', 'PF3U', 'PF5', 'PF3L', 'PF2L', 'PF1AL'};

vvgroup = [1  1  2  2  3  4  5  5  5  6  6  6  6  6  6  6  6  6  6  6 ...
     6  6  6  6  7  7  7  7  7  8  8  8  8  9  9  9  9  9  9  9  9  9 ...
     9  9  9  9  9  9  9  9  9  9  9  9  9  9  9 10 10 10 11 11 11 11 ...
    11 12 12 13 13 13 14 15 16 17 18 18 18 19 19 20 20 20 20 20 21 21 ...
    21 22 22 22 22 22 22 22 22 22 22 22 22 22 22 22 22 22 22 22 22 22 ...
    22 22 22 23 23 23 23 24 24 24 24 24 25 25 25 25 25 25 25 25 25 25 ...
    25 25 25 25 25 26 26 26 27 28 29 29 30 30 31 31 32 32 33 34 35 36 ...
    37 38 39 40];  

% parallel connections, current fraction determined by vessel resistance
vvcirc = 1:max(vvgroup);
resv = tok_data_struct.resv;
Pvv = zeros(length(vvcirc), max(vvcirc)); 
for ii=1:max(vvcirc)
  idx1=find(vvgroup==ii);
  sum_rinv = sum(1./resv(idx1));
  Pvv(idx1,ii)= 1./resv(idx1)/sum_rinv;  
end

ncx = max(cccirc);
nvx = max(vvcirc);
iicx = 1:ncx;
iivx = ncx + (1:nvx);

circ = variables2struct(ccnames, cccirc, vvcirc, vvgroup, Pcc, Pvv, iicx, iivx);
