% This script runs and times several example quadratic programs using the
% Matlab built-in quadprog.m and direct-substitution quadratic program,
% quadprogds.m

% small QP: built-in is faster (0.0008s vs. 0.0019)
run_and_time_qps(0, 10, 2, 3)      

% medium QP: comparable speed (0.21s vs 0.21s)
run_and_time_qps(0, 1000, 10, 30)  

% large QP: direct-substitution is faster (1.6s vs 0.9s)
run_and_time_qps(0, 2000, 400, 30) 

function run_and_time_qps(seed, n, neq, nineq)
% seed = rng seed (for reproducibility)
% n = number of optimization variables
% neq = number of equality constraints
% nineq = number of inequality constraints 
rng(seed);

sqrtH = rand(n,n);
H = sqrtH' * sqrtH;
f = rand(n,1); 
A = rand(nineq,n);
b = rand(nineq,1);
Aeq = rand(neq,n);
beq = rand(neq,1);
lb = [];
ub = []; 
x0 = []; 
options = optimoptions('quadprog', 'display', 'off');

% Run standard quadprog
tic
x1 = quadprog(H,f,A,b,Aeq,beq,lb,ub,x0,options);
t1 = toc;

% Run quadprog with direct substitution
tic
x2 = quadprogds(H,f,A,b,Aeq,beq,lb,ub,x0,options);
t2 = toc; 

% printout
assert(norm(x1-x2) < sqrt(eps)*norm(x1), 'Solver values are different')
fprintf('\n%d vars with %d equality constraints and %d inequality constraints:\n', n, neq, nineq)
fprintf('Built-in quadprog took %.4fs\n', t1);
fprintf('Direct-substitution quadprog took %.4fs\n', t2);
fprintf('Norm of difference between the two methods: %e\n', norm(x1-x2))
end
