function [Cext] = Cext2DOpticalTheorem(ki, r_sca, a, l)
% Cext2DOpticalTheorem(ki, r_sca, a)
%   Compute the extinction cross-section by means of the 2D Optical Theorem.
%   Works very much faster than the functions that perform integration.
%   Arguments:
%   ki - incident wave vector;
%   r_sca - radius-vectors of scatterers: column=scatterer number, row=coordinate
%   a - amplitudes of scattered Hankel waves: column=scatterer number, 
%       row=multipole number (see MultipleRodsHankelAmplitudes)
%   l - Column of multipoles to take into consideration.

  [phi_ki, k] = cart2pol(ki(1), ki(2));
  C = repmat(sqrt(2 / pi / k) ...
          .* exp(-1i .* pi/4) ...
          .* exp(1i .* l .* phi_ki) ...
          .* 1i.^(-l), ...
          1, size(a, 2));
    
  % Repeat ki for each scatterer
  K = repmat(ki, 1, size(r_sca, 2));
  
  % dot-product ki with radius-vectors of each scatterer, exponentiate,
  % and repeat the result for each multipole
  E = repmat(exp(-1i .* dot(K, r_sca)), size(a, 1), 1);
  
  Fforward = sum(sum(C .* E .* a));
  Cext = 2 * sqrt(pi / k) .* (imag(Fforward) - real(Fforward));
end
