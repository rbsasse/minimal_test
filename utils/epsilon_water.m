function epsilon = epsilon_water(wavelength, temperature)
% after Measurement of the refractive index of distilled water from the near-infrared region to the ultraviolet region
% Appl Opt. 2007 Jun 20;46(18):3811-20

% plot(lambda, [epsilon_water(lambda, 10), epsilon_water(lambda, 22), epsilon_water(lambda, 30)], lambda, lambda*0+1.33^2, '--')

if nargin < 2
    temperature = 20;
end

coefs = [0.5672526103, 0.1736581125, 0.02121531502, 0.1138493213, 0.005085550461, 0.01814938654, 0.02617260739, 10.73888649;
0.5684027565, 0.1726177391, 0.02086189578, 0.1130748688, 0.005101829712, 0.01821153936, 0.02620722293, 10.69792721;
0.5689093832, 0.1719708856, 0.02062501582, 0.1123965424, 0.005110301794, 0.01825180155, 0.02624158904, 10.67505178;
0.566695982, 0.1731900098, 0.02095951857, 0.1125228406, 0.005084151894, 0.01818488474, 0.02625439472, 10.73842352];

temps = [ 19          20         21.5           24];

Ni = 4;

wavelength = wavelength(:); % ensure column format
epsilon = 0*wavelength;

for(ii=1:Ni) % 4 terms
    Ai = interp1(temps, coefs(:,ii), temperature, 'spline');
    Li = interp1(temps, coefs(:,ii+4), temperature, 'spline');
    epsi = sellmeier(wavelength*1e-3, Ai, Li);
    epsilon = epsilon + epsi;
end

epsilon = 1 + epsilon;
end

function epsilon = sellmeier(wavelength, Ai, Li)
epsilon = Ai * wavelength.^2 ./ (wavelength.^2 - Li);
end

%   coef          19          20         21.5           24
% 1   A1  0.56725261  0.56840276  0.568909383  0.566695982
% 2   A2  0.17365811  0.17261774  0.171970886  0.173190010
% 3   A3  0.02121532  0.02086190  0.020625016  0.020959519
% 4   A4  0.11384932  0.11307487  0.112396542  0.112522841
% 5   L1  0.00508555  0.00510183  0.005110302  0.005084152
% 6   L2  0.01814939  0.01821154  0.018251802  0.018184885
% 7   L3  0.02617261  0.02620722  0.026241589  0.026254395
% 8   L4 10.73888649 10.69792721 10.675051780 10.738423520