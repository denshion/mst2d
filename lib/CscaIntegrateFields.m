function [Csca] = CscaIntegrateFields(ki, r_sca, a, l, Rint, polarization, gg_)
% CscaIntegrateFields(ki, r_sca, a, l, Rint, polarization, [gg = 100])
%   Compute the scattering cross-section by integrating the scattered fields
%   along a circle with radius Rint with the center at the origin.
%   Arguments:
%   ki - incident wave vector;
%   r_sca - radius-vectors of scatterers: column=scatterer number, row=coordinate
%   a - amplitudes of scattered Hankel waves: column=scatterer number, 
%       row=multipole number (see MultipleRodsHankelAmplitudes)
%   l - multipoles taken into account
%   Rint - radius to integrate at (respective to the origin). Should be large,
%   so that Rint * norm(ki) >> 1.
%   polarization - direction of the electric field vector in cylindrical
%     coordinates:
%     [0, -1, 0] for TM, [0, 0, 1] for TE;
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
      warning('Rint is too small: some scatterers lie out of the circle.');
    end
  end

  k = norm(ki);
  
  % kR << N
  N = fix(gg .* k .* Rint);
  angles = linspace(0, 2 * pi, N);
  integrand = zeros(1, numel(angles));
  
  for iter = 1:numel(angles)
    phi = angles(iter);
    [x, y] = pol2cart(phi, Rint);
    r = [x; y];
    r_3d = [x; y; 0];
    % r - r_sca
    rmr_sca = repmat(r, 1, size(r_sca, 2)) - r_sca;
    [Phi, R] = cart2pol(rmr_sca(1, :), rmr_sca(2, :));
    Hankels = bsxfun(@besselh, l, k .* R);
    ComplexExponents = exp(1i .* l * Phi);
    ScatteredField = sum(sum(a .* Hankels .* ComplexExponents));
    
%    % Wave vector is ki, in cylindrical coordinates it is:
%    ki_pol = [cos(phi), sin(phi); -sin(phi), cos(phi)] * ki;
%    Hi_direction = cross([ki_pol; 0] ./ k, polarization);
    
    % Wave vector is along e_r
    Hs_direction = cross([1, 0, 0], polarization);
  
%    Ei = exp(1i .* dot(ki, r)) .* polarization;
%    Hi = exp(1i .* dot(ki, r)) .* h_direction;
    Es = ScatteredField .* polarization;
    Hs = ScatteredField .* Hs_direction;
    
    integrand(iter) = dot(r_3d, real(cross(Es, conj(Hs))));
  end
  
  Csca = trapz(angles, integrand);
end