%  inputs: settings struct
%          if settings.plotlevel >= 2, makes a plot of the plasma scalars
% 
%  outputs: the plasma_scalars struct that has fields:
%     ip, li, wmhd, Rp each with subfields (Time, Data, Units). 
%     Units is only used for plotting, i.e. don't change units 
%     from A to kA and expect the code to intelligently adapt. 
%
function plasma_params = define_plasma_params(settings)


pp.Ip.Time = [0 0.05 0.15 0.25 0.4 1]';
pp.Ip.Data = [0 2.2 4.7 5.8 7.7 7.7]' * 1e5;

pp.Wk.Time = [0 0.13 0.27 0.49 0.71 1]';
pp.Wk.Data = [0.01 1.2 3 8 8 3.5]' * 1e4;

pp.Rp.Time = [0 0.03 0.075 0.12 0.15 0.24 0.39 0.55 1]';
pp.Rp.Data = [17 12 8 4.7 3.6 2.5 1.8 1.2 1.2]' * 1e-6; 


npsin = 101;
psin = linspace(0,1,npsin);

pp.pprime.Time = [0 1]';
pp.pprime.Data = [1 1]' * -4 * (-(psin-0.5).^2 + 0.25);  

pp.ttprime.Time = [0 1]';
pp.ttprime.Data = [1 1]' * (1-psin).^1.7;

% plotting
if settings.plotlevel >= 2
  plot_structts(pp, {'Ip','Wk','Rp'}, 2);
  sgtitle('Plasma scalar targets'); drawnow
end 

plasma_params = pp;