function ysm = nonuniform_movmean(x,y,xwidth)
% ========================================================================
% EXAMPLE:
% n = 3;
% x = linspace(-4*pi,4*pi,500)';
% x(40:90) = [];
% rng default  
% y = sin(x)*(1:n) + 0.25*rand(length(x),n);
% ysm = nonuniform_movmean(x, y, 0);
% 
% clf
% plot(x,y, 'o')
% hold on
% set(gca, 'ColorOrderIndex', 1)
% plot(x,ysm,'-', 'linewidth', 2)
% legend('Input Data','Filtered Data')
% shg
% ========================================================================

  n = length(x);
  ysm = nan(size(y));
  
  for i = 1:n
    k = x <= x(i) + xwidth/2 & ...
        x >= x(i) - xwidth/2;
    ysm(i,:) = mean(y(k,:), 1);
  end
end
