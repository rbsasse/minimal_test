function epsAu = epsilon_Au(lambda)
%EPSILON_AU Gold dielectric function
%
% Uses the analytical expression given in Eq. (E.1).
% The exp(-i omega t) convention is assumed.
%
% PARAMETERS:
% - lambda: scalar, vector, or matrix
%           wavelengths in NANOMETERS (nm)
%
% RETURNS: epsilon(lambda) as a complex column vector
%
% DEPENDS: none
%
% FAMILY: user_level, epsilon, utility
%

lambda = lambda(:); % ensure column
eps_infty = 1.54;
lambda_p = 177.5;
mu_p = 14500.0;
A1=1.27;
lambda1=470.0;
mu_p1=1900.0;
A2=1.1;
lambda2=325.0;
mu_p2=1060.0;
phi=-pi/4;

epsAu = eps_infty * (1 - 1 ./ (lambda_p.^2 *( (1./lambda).^2 + 1i./(mu_p.*lambda))))...
    + A1 / lambda1 *(exp(1i*phi)./(1/lambda1-1./lambda-1i/mu_p1)+exp(-1i*phi)./(1/lambda1+1./lambda+1i/mu_p1))...
    + A2 / lambda2 *(exp(1i*phi)./(1/lambda2-1./lambda-1i/mu_p2)+exp(-1i*phi)./(1/lambda2+1./lambda+1i/mu_p2));
end
