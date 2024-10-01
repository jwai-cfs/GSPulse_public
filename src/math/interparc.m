function [x2,y2,s] = interparc(x,y,n,forceloop,mergeit)
% =========================================================================
% Description: 
%  interpolate a curve defined by the points (x,y), such that the returned
%  curve has points equally-spaced by arc length
% 
% Inputs: 
% (x,y) - the input curve coordinates
% n - number of points in the returned curve
% mergeit: include original x,y in x2,y2
% forceloop: force curve to be a loop (ie x(1),y(1) = x(end),y(end)
% 
% Outputs: 
% (x2,y2) - equally-spaced curve
% s - cumulative arc length along the curve
% =========================================================================

if ~exist('mergeit','var'), mergeit = 0; end

x = x(:);
y = y(:);

if forceloop
  if x(1) ~= x(end) || y(1) ~= y(end)
    x(end+1) = x(1);
    y(end+1) = y(1);
  end
end

arclens = cumsum([0; sqrt(diff(x).^2 + diff(y).^2)]);

% evenly distributed arc lengths
s = linspace(0,arclens(end),n+1);

for i = 1:n
    k = find(arclens > s(i), 1) - 1;
    dk = (s(i)-arclens(k)) / (arclens(k+1)-arclens(k)); % remainder
    
    x2(i) = x(k) + dk * (x(k+1) - x(k));
    y2(i) = y(k) + dk * (y(k+1) - y(k));
end

s(end) = [];

% Also include original x,y in outputs x2,y2
if mergeit
    
    v = [x2(:) y2(:) s(:); x(:) y(:) arclens(:)];
    v2 = sortrows(v,3);
    
    s = v2(:,3);            
    [s,iuniq] = unique(s);
    
    x2 = v2(iuniq,1);
    y2 = v2(iuniq,2);    
end
x2 = x2(:);
y2 = y2(:);
s = s(:);  