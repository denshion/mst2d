function [Cext] = CextIntegrateScatteringAmplitude(ki, r_sca, a, l, Rint, gg_)
% CextIntegrateScatteringAmplitude(ki, r_sca, a, l, Rint, [gg = 100])
%   Compute the extinction cross-section by integrating the scattering
%   amplitude along a circle with radius Rint with the center at the origin.
%   Arguments:
%   ki - incident wave vector;
%   r_sca - radius-vectors of scatterers: column=scatterer number, row=coordinate
%   a - amplitudes of scattered Hankel waves: column=scatterer number, 
%       row=multipole number (see MultipleRodsHankelAmplitudes)
%   l - multipoles taken into account
%   Rint - radius to integrate at (respective to the origin). Should be large,
%   so that Rint * norm(ki) >> 1.
%   gg - gg * norm(ki) * Rint points will be taken for integration. This value
%   should be large. Default value: s = 100.

  if nargin > 6
    gg = gg_;
  else
    gg = 100;
  end

  if ~iscolumn(ki)
    error('Incident plane wave wave vector should be a column-vector');
  end
  
  if ~iscolumn(l)
    error('Multipoles should be given as a column-vector');
  end
  
  if size(r_sca, 2) ~= size(a, 2)
    error('Number of scatterers in r_sca and a is not equal');
  end
  
  if size(a, 1) ~= numel(l)
    error('Number of multipoles given ~= number of multipoles selected');
  end
  
  if size(ki, 1) ~= size(r_sca, 1)
    error(['Incident plane wave vector and radius-vectors of scatterers ' ...
    'should have be of the same dimension']);
  end

  for i = 1:size(r_sca, 2)
    if norm(r_sca(:, i)) >= Rint
      warning(sprintf( ...
      'Rint=%e too small: some scatterers lie out of the circle.', Rint));
    end
  end

  k = norm(ki);
  
  if k .* Rint <= gg
    warning(sprintf( ...
    ['Rint=%e too small: norm(ki) * Rint=%e is not much greater than unity, ' ...
    'integral might converge poorly'], Rint, k .* Rint));
  end
  
  % 1 << kR << N
  N = fix(gg .* k .* Rint);
  angles = linspace(0, 2 * pi, N);
  integrand = zeros(1, numel(angles));

  for iter = 1:numel(angles)
    phi = angles(iter);
    
    C = repmat(sqrt(2 / pi / k) ...
            .* exp(-1i .* pi/4) ...
            .* 1i.^(-l) ...
            .* exp(1i .* l .* phi), ...
            1, size(a, 2));
    
    % Generate a k from |k| and angle, and repeat it for each scatterer
    [kx, ky] = pol2cart(phi, k);
    K = repmat([kx; ky], 1, size(r_sca, 2));
    
    % dot-product k with radius-vectors of each scatterer, exponentiate,
    % and repeat the result for each multipole
    E = repmat(exp(-1i .* dot(K, r_sca)), size(a, 1), 1);

    F = sum(sum(C .* E .* a));
    
    integrand(iter) = -sqrt(Rint) .* real(...
       exp(1i .* k .* Rint .* (1 - cos(phi))) .* ...
       (1 + cos(phi)) .* ...
       F);
  end
  
  Cext = trapz(angles, integrand);
end
