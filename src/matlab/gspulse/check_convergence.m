function [converged, convergence_msg, prev_ic, nconv] = check_convergence(prev_ic, kiter, settings, nint, stage,nconv)

converged = false;
convergence_msg = '';

if settings.tol_ic_diff > 0  
  
  % read ic
  for kint = 1:nint
    ic{kint} = stage{kint}.mpcsoln.ic.Data;
  end

  % check converged
  if kiter > 1
    
    meets_ic_check = true;
    
    % check if coil current condition is met (coil current difference
    % between iterations is below tolerance)
    max_ic_rel_diff = 0;
    for kint = 1:nint      
      ic_diff = vecnorm(ic{kint} - prev_ic{kint}, 2, 1)';
      ic_diff = ic_diff / sqrt(length(ic{kint})); 
      ic_scale = settings.ic_max - settings.ic_min;
      max_ic_rel_diff = norm(ic_diff ./ ic_scale);
      meets_ic_check = meets_ic_check && (max_ic_rel_diff < settings.tol_ic_diff);
    end

    % check overall convergence 
    % must meet ic_check for 3 successive iterations to be considered converged
    if meets_ic_check
      nconv = nconv+1;
    else
      nconv = 0;
    end
    converged = nconv >= 3; 
           
    % Format convergence message
    if ~meets_ic_check
      convergence_msg = sprintf("Coil currents diff is above prescribed tolerance: %.2e > %.2e", max_ic_rel_diff, settings.tol_ic_diff);
    else
      convergence_msg = sprintf("Coil currents diff is below prescribed tolerance: %.2e < %.2e, \nCondition has been satisfied for %d/3 successive iterations", max_ic_rel_diff, settings.tol_ic_diff, nconv);
    end
    
    if converged
      convergence_msg = strcat(convergence_msg, '\nMet convergence criterion. Stopping iterations.\n');
    else
      if kiter < settings.niter
        convergence_msg = strcat(convergence_msg, '\nContinuing to next iteration.\n');
      else
        convergence_msg = strcat(convergence_msg, sprintf('\nStopping. Reached max number of iterations %d without converging.\n', settings.niter));
      end
    end
  end
  
  % write ic
  for kint = 1:nint
    prev_ic{kint} = ic{kint};
  end
end
