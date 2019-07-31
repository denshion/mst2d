function [s, t] = lorenz_mie_2d(polarization, ...
                                n, ...
                                l, ...
                                x, ...
                                admittance)
% [s, t] = lorenz_mie_2d(polarization, n, l, x, [admittance])
% Calculate the Lorenz-Mie coefficients of a dielectric rod
%   polarization - 'TE' or 'TM';
%   n - rod refraction index divided by filling medium
%       refraction index (note that it's equal to
%       relative_admittance for non-magnetic media)
%   l - multipole number;
%   x - size parameter kR (wave vector multiplied by radius);
%   admittance - rod admittance divided by filling medium admittance;\
% The function is GPU-optimized.
% TODO: allow negative x.

  if nargin < 5
    admittance = n;
  end
  if strcmp(polarization, 'TE')
    p_e = -1;
    p_i = -1 ./ admittance;
  elseif strcmp(polarization, 'TM')
    p_e = 1;
    p_i = admittance;
  else
    error('polarization should be TE or TM')
  end

  nx = x .* n;
  % GPU-optimized. GPU only allows non-negative NU, so we process the
  % negative case explicitly here.
  y0 = sign(l) .^ (l) .* bsxfun(@bessely, abs(l), x);
  y1 = sign(l+1) .^ (l+1) .* bsxfun(@bessely, abs(l + 1), x);
  j0 = sign(l) .^ (l) .* bsxfun(@besselj, abs(l), x);
  j1 = sign(l+1) .^ (l+1) .* bsxfun(@besselj, abs(l + 1), x);
  jn0 = sign(l) .^ (l) .* bsxfun(@besselj, abs(l), nx);
  jn1 = sign(l+1) .^ (l+1) .* bsxfun(@besselj, abs(l + 1), nx);

  alpha = (p_e .* jn0 .* (l.*y0 ./ x -  y1) - p_i .* y0  .* (l.*jn0./nx - jn1)) ...
       ./ (p_i .* j0  .* (l.*jn0./nx - jn1) - p_e .* jn0 .* (l.*j0 ./ x -  j1));
  
  s = 1 ./ (-1 + 1i .* alpha);
  if nargout > 1
    t = 1i .* s .* ...
       (p_e .* j0 .* (l.* y0./ x -  y1) - p_e .*  y0 .* (l.*j0./x - j1)) ...
    ./ (p_i .* j0 .* (l.*jn0./nx - jn1) - p_e .* jn0 .* (l.*j0./x - j1));
  end
end