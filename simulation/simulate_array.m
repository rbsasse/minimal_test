function  xsec = simulate_array(rho_dye, N, tensor_type, angle, wavelength, medium)
%simulate_array Simulates the optical response of a square array cluster
%
% High-level function to wrap a full simulation 
% (to be called by a script with varying parameters)
%
% PARAMETERS:
% - rho_dye: surface density of dipoles
% - N: number of dipoles along one side
% - tensor_type: type of tensor [1 'isotropic', 2 'radial', 3 'random-flat', 4 'flat isotropic']
% - angle: angle of incidence
% - polarisation: polarisation
% - wavelength: vector of wavelengths
% - medium: refractive index of surrounding medium
%
% RETURNS: matrix of cross-section spectra
%
% DEPENDS: alpha_bare, cluster_array, spectrum_dispersion
%
% FAMILY: user_level, simulation
%

pitch = 1/sqrt(rho_dye);

switch tensor_type
    case 1 % isotropic
        a = 1; b=1; c=1;
    case 2 % normal
        a = 0; b=0; c=1;
    case 3 % flat-isotropic
        a = 1; b=1; c=0;
    otherwise
        warning('Unexpected tensor type.')
end

alpha_inf=9.6e-39;
alpha_k = [5.76e-38];
lambda_k = [526];
mu_k = [10000];

alphabar = alpha_bare(wavelength, alpha_inf, alpha_k, lambda_k, mu_k);
cl = cluster_array(N, pitch, a, b, c);

material.wavelength = wavelength;
material.alpha = alphabar;
% cater for dispersive media
if isnumeric(medium) && length(medium) == 1
    material.medium = material.wavelength*0 + medium;
else % all(ishandle(medium)) || all(ischar(medium))
    material.medium = feval(medium, wavelength);
end

Incidence = [0;angle;0];
xsec = spectrum_dispersion(cl, material, Incidence);
end
