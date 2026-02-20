function Sm = get_spline_transform(t, settings, nu, nx)

% Spline transformation
Ic_scale = 10e3;
Vc_scale = 1e3;

% compute move-blocking (spline basis) matrices
% for use as in eqns A.42 and A.43
if settings.use_spline_basis
  n_nospline = 2; % don't spline compress the first 2 time steps
  t_spline = t(n_nospline+1:end); % time points for which to compress
  n_basis_pts = round(numel(t_spline) * settings.spline_basis_ratio);
  t_basis = linspace(t_spline(1), t_spline(end), n_basis_pts);  
  [~,bm] = spline_basis(t_basis, nan(n_basis_pts,1), t_spline);  
  bm = blkdiag(eye(n_nospline), bm); % basis for one time trace
  Sm_u = kron(bm, eye(nu)); % so that [u0;u1,...] = Sm_u*p
else
  Sm_u = eye(nu*numel(t));
end

Sm_x = eye(nx);

% build final Sm including scaling so that [x1;u0;u1..] = Sm*x
Sm = blkdiag(Ic_scale * Sm_x, Vc_scale*Sm_u);
Sm = sparse(Sm);