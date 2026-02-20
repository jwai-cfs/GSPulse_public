function mat = matrix_diff_first_derivative(times)

N = length(times);
mat = zeros(N-1,N);

for i = 2:N
  dt = times(i) - times(i-1);
  mat(i-1,i-1) = 1 / dt;
  mat(i-1,i)   = -1 / dt;
end
