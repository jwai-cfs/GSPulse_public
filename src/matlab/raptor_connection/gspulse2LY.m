function LYt = gspulse2LY(soln, L)
% Convert GSPulse outputs to Raptor inputs
%
% This function takes the solution from a GSPulse run and transforms the
% output equilibria to the LYt format, which is used as inputs for RAPTOR.

eqs = soln.eqs;
t = soln.t;
N = length(t);

% sometimes eqs have additional fields, strip the additional fields so that
% they are all consistent (otherwise LYt assignment fails)
fds = fieldnames(eqs{1});

for i = 1:N
  fds2remove = setdiff(fieldnames(eqs{i}), fds);
  eqs{i} = rmfields(eqs{i}, fds2remove);
end

% copy to LYt struct
for i = 1:N
  eqs{i}.rx = L.G.rx;
  eqs{i}.zx = L.G.zx;
  eqs{i}.rg = L.G.rx;
  eqs{i}.zg = L.G.zx;
  eqs{i}.t = t(i);
end
LYt = cell2mat(eqs);