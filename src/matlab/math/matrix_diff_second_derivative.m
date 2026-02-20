function mat = matrix_diff_second_derivative(times)

N = length(times);
mat = zeros(N-2, N);

for k = 3:N
  dt_k_km1 = times(k) - times(k-1);
  dt_km1_km2 = times(k-1) - times(k-2);
  dt = (dt_k_km1 + dt_km1_km2)/2;

  mat(k-2,k)   = 1 / (dt * dt_k_km1);
  mat(k-2,k-1) = -1 / (dt*dt_k_km1) - 1 / (dt*dt_km1_km2);
  mat(k-2,k-2) = 1 / (dt * dt_km1_km2);
end
