function L = load_anamak_L_fbt(P)


% These are variations from the default "fbtana.m" settings that are not
% captured by the GUI inputs already: 
P.gsxe = 3;       % force Mxx calculation
P.selu = 'e';     % use eigenmodes for vessel currents
P.ifield = true;  % compute grid (Br,Bz)

% Call FBT to generate L
params = [fieldnames(P) struct2cell(P)]';
L = fbt('ana', 0, 0, params{:});
