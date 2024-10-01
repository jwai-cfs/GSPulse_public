function mp = unwrap_mpp(mpp, nz, nr)
% =========================================================================
% Description: 
%  TokSys often stores the mpp (plasma grid to plasma grid mutuals) in a
%  compressed format that has dimension [nz*nr x nr]. This unwraps it to
%  the "full" format that has dimension [nz*nr x nz*nr]. 
% 
% Inputs: 
%   mpp - plasma mutuals compressed format
%   nz - number of vertical grid points
%   nr - number of radial grid points
% 
% Outputs: 
%   mp - plasma mutuals, uncompressed format
%
% =========================================================================

ngrids = size(mpp,1);
[ridx, zidx] = meshgrid(1:nr, 1:nz);
mp = zeros(ngrids);

for j=1:ngrids
  for k = 1:ngrids
    n = mod(k-1,nz)+1 - zidx(j);
    m = floor((k-1)/nz)+1;
    mp(k,j) = mpp((m-1)*nz+abs(n)+1,ridx(j));
  end
end