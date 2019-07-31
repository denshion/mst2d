function [a] = RodsHankelAmplitudes(k, n, rods_radii, rods_positions, l, excitation, polarization)
% a = RodsHankelAmplitudes(k, n, rods_radii, rods_positions, l, polarization)
%  Compute the Hankel amplitudes after plane wave scattering on multiple rods
%  Arguments:
%  k - incident wave vector
%  n - refraction indices of the rods
%  rods_radii - radii of the rods
%  rods_positions - positions of the rods
%  l - Column of multipoles to take into consideration
%  polarization - 'TE' or 'TM'
% The amplitudes are given in a matrix. The rows correspond to multipoles
% as given in the l column; the columns correspond to numbers of the rods.
  
  rods_count = numel(rods_radii);
  M = numel(l);
  
  if size(rods_positions, 2) ~= rods_count
    error('numel(positions) != numel(radii)')
  end
  
  if numel(n) ~= rods_count
    error('numel(n) != numel(radii)')
  end
  
  if ~iscolumn(l)
    error('Multipoles should be given as a column-vector');
  end
  
  lmm = bsxfun(@minus, l', l);
          
  A = eye(M .* rods_count, M .* rods_count);
  u = zeros(M .* rods_count, 1);
  for i = 1:rods_count    
    s = lorenz_mie_2d(polarization, n(i), l, k .* rods_radii(i));
    
    % Plane wave excitation of a rod at r(i)
    block_i = ((i - 1) * M + 1):(i * M);
    u(block_i) = excitation(:, i) .* s;

    % Rods cross-interactions:
    for j = 1:rods_count
      block_j = ((j - 1) * M + 1):(j * M);
      if i ~= j
        r_ij = rods_positions(:, i) - rods_positions(:, j);
        [phi_ij, R_ij] = cart2pol(r_ij(1), r_ij(2));
        L_ij = besselh(lmm, k .* R_ij) .* exp(-1i .* lmm .* phi_ij);
        A(block_i, block_j) = -diag(s) * L_ij;
      end
    end
  end
  
  a_vector = A \ u;
  
  % Unwrap Hankel Amplitudes vector to matrix:
  % multipoles in columns, rod numbers in rows
  a = reshape(a_vector, M, rods_count);
end