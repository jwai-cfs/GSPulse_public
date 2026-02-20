function out = get_inequality_constraint(...
  settings, ny, nu, nx, Nlook, M, cv, kinterval)

% add inequality constraints for the current limits
ymin = -inf(ny,1);
ymax = inf(ny,1);
ymin(cv.iy.ic) = settings.ic_min;
ymax(cv.iy.ic) = settings.ic_max;
yminhat = repmat(ymin, Nlook, 1);
ymaxhat = repmat(ymax, Nlook, 1);  

if kinterval > 1   
  % except for the first stage, constraints on the first 2 timesteps are
  % enforced separately (implicitly via e1/e2/x1/u0/u1), so we need to
  % remove the constraints from the first 2 timesteps here or 
  % else the QP can be overconstrained and not pass numerical checks. 
  yminhat(1:2*ny) = -inf;
  ymaxhat(1:2*ny) = inf; 
end

Aineq = [-M; M];  % eqn A.31
% bineq must be updated at each iteration, since d changes across
% iterations. For now, specify nan. Note the equation is: 
% bineq = [ymaxhat + d - rhat; -yminhat - d + rhat];
bineq = nan(size(Aineq,1), 1);

% add inequality constraints for the voltage limits
vmax = repmat(settings.vmax, Nlook, 1);
vmin = repmat(settings.vmin, Nlook, 1);

if kinterval > 1   
  % Except for the first stage, u0 and u1 are constrained already.
  % Can fail numerical tolerance checks if inequality constraints also 
  % specified here.
  vmax(1:nu*2) = inf(nu*2,1);
  vmin(1:nu*2) = -inf(nu*2,1);
end

I = [zeros(Nlook*nu, nx) eye(Nlook*nu)];
Aineq = [Aineq; I; -I];
bineq = [bineq; vmax; -vmin];

out = variables2struct(Aineq, bineq, ymaxhat, yminhat);
