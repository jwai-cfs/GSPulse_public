function out = get_equality_constraint(...
  cv, optimization_signals,interval_time_data,init,nu,nx,ny,Nlook,M,kinterval)

% put current and voltage constraints into the form Aeq*phat = beq 
[~,~,~,targs,islocked] = assemble_weights_and_targets(optimization_signals);

A = {}; % will hold equality constraint matrices for each source for Aeq
b = {}; 

% constraint currents
% --------------------

% form coil current constraints from waveform, applies to t=2...N
ilock = logical(structts2vec(islocked, cv.y_names, interval_time_data.t1N));
ilock(1:ny) = 0; % remove constraint at t=1
A{1} = M(ilock,:); % equality constraint version of eqn A.34, with e=0, is M(ilock,:) * phat = -d(ilock)
b{1} = nan(sum(ilock),1); % b := -d(ilock) cannot be precomputed since d is update each iteration
idx_current_locks_t2N = ilock;


% constrain x1 from initial condition (x1 := coil AND vessel currents at t=1) 
if kinterval == 1  
  % init is the global initial condition, which is only a source of
  % constraints for the very first interval. Parse init and determine which
  % vals are free vs locked, nans indicate free parameters. 
  vals = init.x1;
  ilock = ~isnan(vals);
  tmp_ = eye(nx,nx+nu*Nlook);
  A{2} = tmp_(ilock,:); 
  b{2} = vals(ilock);
else
  % for subsequent intervals, the initial condition constraint is
  % determined by the solution of the previous interval, in order to have
  % continuity of overlapping stages. See Appendix A.4. Because it depends
  % on the solution, we cannot assign values for b here and these must be
  % updated per iteration. 
  A{2} = eye(nx, nx+nu*Nlook);
  b{2} = nan(nx,1);   % placeholder to update
end


% constrain voltages
% -------------------
% parse constraint from waveform specified constraints
ilock = logical(structts2vec(islocked, {'voltage'}, interval_time_data.t0Nm1));
vals = structts2vec(targs, {'voltage'}, interval_time_data.t0Nm1);

% Now also add u0 and u1 constraint (voltage at t=0, t=1) from initial 
% condition. We overrwrite the waveform vals here, but this is safe b/c
% consistency for the 2 sources was enforced in validate_gspulse_inputs.m
if kinterval == 1
  % first interval, constraints come from global initial condition
  if Nlook == 1  % special case, constrain only u0
    ilock(1:nu) = ~isnan(init.u0);
    vals(1:nu) = init.u0;
  else            % standard case, constrain u0 and u1
    ilock(1:(2*nu)) = ~isnan([init.u0; init.u1]);
    vals(1:(2*nu))  = [init.u0; init.u1];
  end
else
  % for subsequent intervals, the initial condition is updated each
  % iteration, to enforce continuity of overlapping stages. Assign
  % placeholder values to be updated each iteration. 
  if Nlook == 1
    ilock(1:nu) = true; 
    vals(1:nu) = nan(nu,1); 
  else
    ilock(1:(2*nu)) = true;
    vals(1:(2*nu)) = nan(2*nu,1);
  end
end
tmp_ = [zeros(nu*Nlook,nx) eye(nu*Nlook)];
A{3} = tmp_(ilock,:);
b{3} = vals(ilock);

% write outputs
% --------------
Aeq = vertcat(A{:});
out = variables2struct(Aeq, A, b, idx_current_locks_t2N);
