function soln = GSPulse(gsp_inputs)
% =========================================================================
% Description: 
%  main body of the GSPulse algorithm, and called by run_pulse.m
%
% Inputs: 
%  gsp_inputs - inputs specified by a users config files. For info on the
%               config files, try help run_pulse.m
% Outputs: 
%  soln - GSPulse output solution struct. Contains:
%    t - timebase for equilibria
%    eqs - cellarray of equilibria
%    signals - waveforms for all of the signals computed by the optimizer
%    stage - individual stage/interval solutions for the equilibria, if
%            using multiple intervals to solve for the optimization
%    shapes - shapes struct, see help _define_shapes.m
%    plasma_params - plasma_params struct, see help _define_plasma_params.m
%    targs - targs struct, possibly with targs.psibry updated from the
%    Ejima equation. See help _define_targs.m
%
% Additional info: 
% A written description of the GSPulse algorithm is available in
% docs/gspulse_algorithm.pdf
%
% =========================================================================

% read inputs
settings = gsp_inputs.settings; 
L = gsp_inputs.L; 
tok = gsp_inputs.tok;
shapes = gsp_inputs.shapes;
plasma_params = gsp_inputs.plasma_params;
targs = gsp_inputs.targs;
weights = gsp_inputs.weights;
init = gsp_inputs.init;
model_offset = gsp_inputs.model_offset;

% read timing parameters
settings.timedata = extract_timedata(settings.interval_t);

% dimensions check: transpose data if necessary
weights.wts = check_structts_dims(weights.wts);
weights.dwts = check_structts_dims(weights.dwts);
weights.d2wts = check_structts_dims(weights.d2wts);
targs = check_structts_dims(targs);
plasma_params = check_structts_dims(plasma_params);
shapes = check_structts_dims(shapes);
model_offset = check_structts_dims(model_offset);

nint = settings.timedata.n_intervals;  % number intervals
stage = cell(nint, 1);           % place holder for all the stage solutions

% initialize the first stage
for i = 1:nint
  stage{i}.t      = settings.timedata.interval(i).t;
  stage{i}.pcurrt = initialize_pcurrt(tok, shapes, plasma_params, stage{i}.t);
  stage{i}.eqs    = initialize_eqs(i, settings, shapes, stage{i}.pcurrt, tok, init);
  stage{i}.ws     = []; 
end
stage{1}.init = init;

% compute target psibry from Ejima equation
if strcmpi(settings.specify_psibry_mode, 'ejima')
  [t, eqs] = merge_stage_solutions(stage, settings.timedata);
  targs.psibry = compute_target_psibry(eqs, t, shapes, plasma_params, tok);         
end

% pre-compute as much of the MPC optimization problem as possible, since
% some of the big matrices only needed to be computed once and don't change
% across GS iterations
for i = 1:nint
  fprintf('Precomputing MPC matrices (interval %d/%d) ...\n', i, nint)
  stage{i}.config = mpc_config(i, tok, shapes, settings, weights, targs);
end

% begin main solver loop - loop across Grad-Shafranov convergence iterations
for kiter = 1:settings.niter
  fprintf('GS iteration %d of %d\n', kiter, settings.niter)

  % loop through time intervals, if solving the optimization in multiple
  % intervals
  for kint = 1:nint

    % overlap initial condition of the next stage with the solution from 
    % the previous stage
    if kint > 1
      cv = define_data_indices(settings, targs, tok);
      timedat = settings.timedata.interval(kint);
      stage{kint}.init = define_init_from_prev_stage(timedat, stage{kint-1}, cv, tok);
      stage{kint}.eqs{1} = copyfields(stage{kint}.eqs{1}, stage{kint}.init.eq1, [], 1);
      stage{kint}.eqs{2} = copyfields(stage{kint}.eqs{2}, stage{kint}.init.eq2, [], 1);     
      stage{kint}.pcurrt.Data(1,:) = stage{kint}.eqs{1}.pcurrt(:);
      stage{kint}.pcurrt.Data(2,:) = stage{kint}.eqs{2}.pcurrt(:);
    end
    
    % optimize current evolution (see Appendix B, Step 4: optimize
    % conductor evolution in gspulse_algorithm.pdf)
    [stage{kint}.mpcsoln, stage{kint}.ws] = mpc_update_psiapp(...
      kint, stage{kint}.pcurrt, stage{kint}.config, tok, shapes, ...
      stage{kint}.init, settings, targs, model_offset, stage{kint}.ws, ...
      stage{kint}.eqs);
  end
  
  % plasma Picard iteration (see Appendix B, Step 5: Grad-Shafranov picard 
  % iteration in gspulse_algorithm.pdf)
  for kint = 1:nint
    switch settings.picard_algo
      case 'gspulse'
        [stage{kint}.eqs, stage{kint}.pcurrt] = gs_update_psipla(kint, ...
          stage{kint}.mpcsoln, stage{kint}.pcurrt, tok, plasma_params, settings);
      case 'fbt'
        [stage{kint}.eqs, stage{kint}.pcurrt] = fbt_update_psipla(kiter, ...
          kint, stage{kint}.mpcsoln, stage{kint}.eqs, stage{kint}.pcurrt,...
          tok, plasma_params, settings, L);
    end  
  end

  % update target psibry from Ejima equation
  if strcmpi(settings.specify_psibry_mode, 'ejima')
    [t, eqs] = merge_stage_solutions(stage, settings.timedata);
    targs.psibry = compute_target_psibry(eqs, t, shapes, plasma_params, tok);         
  end

end

% perform post-processing steps and save outputs
[t, eqs, signals] = merge_stage_solutions(stage, settings.timedata);
[eqs, signals] = gsp_postprocessing(t, eqs, signals, tok, settings);
soln = variables2struct(t, eqs, signals, stage, shapes, plasma_params, targs);