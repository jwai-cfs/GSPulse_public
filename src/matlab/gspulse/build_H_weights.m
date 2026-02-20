function [Hp,He] = build_H_weights(cv, Nlook, nx, nu, ny, interval_time_data, ...
  optimization_signals)

[wts, dwts, d2wts] = assemble_weights_and_targets(optimization_signals);

% build u and e transform matrices: eqns A.18-A.20
% intialize to zero
wu = 0; wdu = 0; wd2u = 0; wde = 0; wd2e = 0;
Su = 0; Sdu = 0; Sd2u = 0; Sde = 0; Sd2e = 0; 

% Read weights for the error term, applies for all values of Nlook
t_we = interval_time_data.t1N;
we   = structts2vec(wts, cv.y_names, t_we);
Se = speye(Nlook*ny); 

% Special case, only solving for a single equilibrium.
if Nlook <= 1

  % the optimization vector is only p = [x1; u0] in this case, and u0 does
  % not actually effect the solution (all the entries of F are zero). So we
  % can safely assign some weights to the voltage term u0 (do not need to
  % read them from user). 
  wu = ones(nu,1);
  Su = eye(Nlook*nu);
  
% When solving for at least 2 equilibria the 0th and 1st-derivative terms apply   
elseif Nlook >= 2

  % voltage: read zeroth and first derivative weights
  t_wu = interval_time_data.t0Nm1;    % times at which to interpolate
  t_wdu = (t_wu(1:end-1) + t_wu(2:end)) / 2;   % use midpoint spacing 
  wu  = structts2vec(wts, cv.u_name, t_wu);
  wdu = structts2vec(dwts, cv.u_name, t_wdu);

  % voltage: zeroth-derivative transform matrix is identity
  Su = speye(Nlook*nu); 

  % voltage: first-derivative transformation matrix
  m = matrix_diff_first_derivative(interval_time_data.t0Nm1);
  Sdu = kron(m, eye(nu));
  Sdu = sparse(Sdu); 

  % outputs: read first-derivative weights
  t_wde = (t_we(1:end-1) + t_we(2:end)) / 2;   
  wde = structts2vec(dwts, cv.y_names, t_wde);
  
  % outputs: build first-derivative transformation matrix
  m = matrix_diff_first_derivative(interval_time_data.t1N);
  Sde = kron(m, eye(ny));
  Sde = sparse(Sde); 

  % When solving for 3 or more equilibria, 2nd derivative terms apply
  if Nlook >= 3
    
    % voltage: second derivative quantities
    t_wd2u = (t_wdu(1:end-1) + t_wdu(2:end)) / 2;
    wd2u = structts2vec(d2wts, cv.u_name, t_wd2u);       
    m = matrix_diff_second_derivative(interval_time_data.t0Nm1);
    Sd2u = kron(m, eye(nu));
    Sd2u = sparse(Sd2u);

    % outputs: second derivate quantities
    t_wd2e = (t_wde(1:end-1) + t_wde(2:end)) / 2;
    wd2e = structts2vec(d2wts, cv.y_names, t_wd2e);
    m = matrix_diff_second_derivative(interval_time_data.t1N);
    Sd2e = kron(m, eye(ny));
    Sd2e = sparse(Sd2e);

  end
end

% build quadprog inputs: eqn A.21
He = Se'*diag(we)*Se + Sde'*diag(wde)*Sde + Sd2e'*diag(wd2e)*Sd2e;
He = sparse(He);

Hu = Su'*diag(wu)*Su + Sdu'*diag(wdu)*Sdu + Sd2u'*diag(wd2u)*Sd2u;
Hu = sparse(Hu);
Hx1 = zeros(nx);
Hp = blkdiag(Hx1, Hu); % eqn A.31, Hp = Tup'*Hu*Tup