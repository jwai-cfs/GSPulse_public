function [X,fval,exitflag,output,lambda] = quadprogds(H,f,A,b,Aeq,beq,lb,ub,x0,options)
% =========================================================================
% Quadratic Program Direct Substitution
%
% Usage: 
%  [X,fval,exitflag,output,lambda] = quadprogds(H,f,A,b,Aeq,beq,lb,ub,x0,options) 
%  [X,fval,exitflag,output,lambda] = quadprogds(H,f,A,b,Aeq,beq,lb,ub,ws) 
% 
% Description: 
%   Solves a quadratic program using direct substitution of the equality
%   constraints as opposed to passing the equality constraints to the
%   QP-solver (which uses the Lagrange multiplier method to satisfy
%   equality constraints). This method works by finding the null space of 
%   the equality constraints and solving a smaller QP for a solution within
%   that null space. 
% 
%   This direct substition method has been observed to be faster in some
%   cases, particularly with large problems or large numbers of equality 
%   constraints. If there are no equality constraints (Aeq=[], beq=[]) 
%   then there is no difference between this and quadprog.m
%
% Method: 
%   Finds the null space of the equality constraints. Uses the
%   substitution:
%   
%   X = xeq + null(Aeq) * Xt
%
%   where X is the original solution candidate, xeq is any particular 
%   solution to Aeq*xeq = beq, and Xt is the (lower-dimensional)
%   transformed solution candidate. This relation is plugged everywhere
%   into the original QP to obtain a lower-dimensional QP formulation. 
%
% Restrictions: 
%   Does not support all the syntax flexibility of the Matlab built-in
%   quadprog. 
%
% Author: 
%   Josiah Wai, 1/13/2026
% =========================================================================
if nargin < 10, options = []; end

if isempty(Aeq)
  % direct passthrough
  [X,fval,exitflag,output,lambda] = quadprog(H,f,A,b,Aeq,beq,lb,ub,x0,options);
else
  NA = null(Aeq);

  % Solve the equality constraint equation
  xeq = Aeq\beq;    
  assert(norm(beq - Aeq*xeq) < norm(beq) * sqrt(eps), ...
    'Equality constraints not satisfied')
  
  % coordinate transformation for H
  if ~issymmetric(H)  % ensure H is symmetric
    if verbosity > -1
      warning(message('optim:quadprog:HessianNotSym'))
    end
    H = (H+H')*0.5;
  end
  Ht = NA' * H * NA; % numerically, this step can remove symmetry
  Ht = (Ht + Ht')/2; % ... so enforce symmetry here

  % coordinate transformations for f,A,Aeq,beq
  ft = NA' * (H'*xeq + f);  
  At = A*NA; 
  bt = b - A*xeq;
  Aeqt = [];
  beqt = [];
  
  % coordinate transformation for x0
  if isa(x0, 'optim.warmstart.QuadprogWarmStart')
    % x0 is actually a warmstart object (lol), do nothing
    x0t = x0;
  elseif isempty(x0)  % do nothing
    x0t = [];
  else  % x0 is an initial guess vector, do transform
    x0t = NA'*(x0-xeq);
  end

  % coordinate transformation for (lb,ub), which have to be handled via 
  % (At,Bt) now  
  if ~isempty(lb)
    At = [At; -NA];
    bt = xeq - lb;    
  end
  if ~isempty(ub)
    At = [At; NA];
    bt = xeq + ub;    
  end
  lbt = [];
  ubt = [];

  % pass to quadprog
  [Xt,fval,exitflag,output,lambda] = quadprog(...
    Ht,ft,At,bt,Aeqt,beqt,lbt,ubt,x0t,options);
  
  % coordinate back-transformation of QP solution
  X = xeq + NA * Xt;  
end
