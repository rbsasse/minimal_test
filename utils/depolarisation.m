function chi = depolarisation(x1, x2, x3)
%DEPOLARISATION Computes depolarisation factors for ellipsoids
%
% General form using numerical integration; scaled internally to help convergence.
% For spheroids, it would make sense to use analytical formulas and vectorise the function.
%
% PARAMETERS:
% - x1: axis length
% - x2: axis length
% - x3: axis length
%
% RETURNS: scalar
%
% DEPENDS: none
%
% FAMILY: low_level, polarizability
%

V = x1*x2*x3;
x2 = x2/x1;
x3 = x3/x1;
integrand = @(q) 1./((1 + q) .* sqrt((q + 1) .* (q + x2^2) .* (q + x3^2)));

I1 = quadgk(integrand, 0, inf);
chi = V/2 * I1 / x1^3;

end
