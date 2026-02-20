function [wsout, ok] = solve_qp(H, f, Aineq, bineq, qpsolver, verbose, tol, ws)

H = (H+H')/2;
if isempty(ws)
  recalculate_ws = 1;
  x_init = -H\f; % unconstrained solution as warmstart guess
else
  recalculate_ws = 0; 
end

switch qpsolver

  % --------------- quadprog solver ------------
  case 'quadprog'
    
    if recalculate_ws
      if verbose
        display='iter';
      else
        display='none';
      end
      qpopts = optimoptions('quadprog', ...
                            'Algorithm', 'active-set', ...
                            'Display', display, ...
                            'OptimalityTolerance', tol);
      ws = optimwarmstart(x_init, qpopts);
    end

    
    % call the solver
    [wsout, ~, exitflag] = quadprog(H, f, Aineq, bineq, [], [], [], [], ws);
    ok = exitflag == 1;
    
  % ---------------- scs solver ---------------
  case 'scs'
    
    data = struct;
    if recalculate_ws
      data.x = x_init;
    else
      data.x = ws.X; data.y = ws.y; data.s = ws.s;
    end
    
    % set up data 
    data.P = H;
    data.c = f;
    data.A = Aineq;
    data.b = full(bineq); 
    cone.l = size(data.A,1); % cone solver auxiliary info
    
    % solver settings
    scsopts = struct('eps_abs', tol, 'eps_rel', tol,'eps_infeas',tol,...
      'verbose',double(verbose>1));

    % solve
    [x, y, s, info] = scs(data, cone, scsopts);
    
    % assign solution and lagrange multipliers to warmstart object
    wsout.X = x;
    wsout.y = y;
    wsout.s = s;
    exitflag = info.status_val;
    ok = info.status_val == 1;

  % ---------------- cvxopt via python ---------------
  case 'python_cvxopt'
   
    % specify warmstart values - this is a placeholder operation, since in 
    % testing with CVXOPT, initialization with a known solution is actually
    % slower (even for general QP problems) TODO: debug this and actually
    % use the warmstart. 
    if recalculate_ws
      x0 = x_init;
    else
      x0 = ws.X;
    end   
    
    % convert any sparse arrays to full
    [H,f,Aineq,bineq,x0] = deal(full(H), full(f), ...
      full(Aineq), full(bineq), full(x0));
    qp = variables2struct(H, f, Aineq, bineq, x0, tol);      

    infile = [tempname() '.mat'];
    outfile = [tempname() '.mat'];   
    pyfile = [getenv('GSROOT') '/src/python/gspulse/qp_solve.py'];
    
    % different call signatures for running python depending on whether
    % environment is matlab or octave
    is_octave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
    if is_octave      
      
      % save inputs to file           
      save('-v7', infile, 'qp')  

      % run python script to solve QP
      [~, exitflag] = python(pyfile, infile, outfile); % call python to solve

    else  % matlab 

      % save inputs to file
      save(infile, 'qp') 

      % set python environment
      uv_python_env = [getenv('GSROOT') '/src/python/.venv/bin/python'];
      pyenv('Version', uv_python_env);

      % run python script to solve QP
      pyrunfile([pyfile ' ' infile ' ' outfile]);
    end
    
    % cleanup
    wsout = load(outfile);
    wsout.X = wsout.X(:);
    delete(infile, outfile);
    ok = 1; 

  % ---------------- unknown solver ---------------
  otherwise
    error('unknown QP solver %s',qpsolver);
end

if ~ok
  error('%s did not converge with exitflag=%d\n', qpsolver, exitflag);
end
