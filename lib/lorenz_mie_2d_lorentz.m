function [x0, y0, A] = lorenz_mie_2d_lorentz(polarization, ...
                                             n_rod, ...
                                             n_filling, ...
                                             relative_admittance ...
                                            )
% [x0, y0, A] = lorenz_mie_2d_lorentz(polarization, ...
%                                     n_rod, ...
%                                     [n_filling, ...
%                                     [relative_admittance]])
%   Approximates the dipole scattering Lorenz-Mie coefficient with a
%   Lorentzian response
%   a0(x) = A ./ ((x - x0) + 1i .* y0).
%   The function returns the resonance frequency x0,
%   the broadening factor y0, and the amplitude A = -1i * y0.
%   polarization - 'TE' or 'TM';
%   n - rod refraction index divided by filling medium
%       refraction index (note that it's equal to
%       relative_admittance for non-magnetic media)
%   n_filling - refraction index of the filling medium (vacuum by default);
%   relative_admittance - rod admittance divided by filling medium admittance;
  if nargin < 3
    n_filling = 1;
  end
  
  if nargin < 4
    relative_admittance = n_rod ./ n_filling;
  end
  
  if strcmp(polarization, 'TE')
    p_e = -1;
    p_i = -1 ./ relative_admittance;
  elseif strcmp(polarization, 'TM')
    p_e = 1;
    p_i = relative_admittance;
  else
    error('polarization should be TE or TM')
  end
  
  n = n_rod ./ n_filling;
  
  tf = @(t) p_e .* besselj(0, n_rod .* t) .* bessely(1, n_filling .* t) ...
          - p_i .* bessely(0, n_filling .* t) .* besselj(1, n_rod .* t);
  
  
  % Quickly sample the function to roughly find all zeros in some "reasonable" region
  if n < 3.5
    estimate_formula = pi ./ 2 ./ (n - 1);
    estimate_min = max(estimate_formula - 0.5 - 0.1, 0);
    estimate_max = estimate_formula + 0.5 + 0.1;
  else
    estimate_min = pi ./ (n - 1) .^ 0.8;
    estimate_max = pi ./ 4 ./ (n - 1) .^ 1.2;
  end
  samplex = linspace(estimate_min, estimate_max, 100);
  
  % Find all points where function changes sign:
  sign_change_points = nonzeros(samplex(find(diff(sign(tf(samplex))))));
  
  % Take the minimal one and find the zero more precisely:
  zero_guess = min(sign_change_points);
  x0 = fzero(tf, zero_guess);

  if nargout > 1
    nx0 = n .* x0;
    v_i = 1 ./ n_rod;
    v_e = 1 ./ n_filling;
    
    y0 = ( ...
        p_i .* besselj(0, x0) .* besselj(1, nx0) ...
      - p_e .* besselj(0, nx0) .* besselj(1, x0) ...
    ) ./ ( ...
        (p_e ./ v_e - p_i ./ v_i) .* bessely(0, x0) .* besselj(0, nx0) ...
      - (p_e ./ v_i - p_i ./ v_e) .* bessely(1, x0) .* besselj(1, nx0) ...
    );
  end
  
  if nargout > 2
    A = -1i .* y0;
  end
end