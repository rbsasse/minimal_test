function  xsec = simulate_shell(rho_dye, tensor_type, position_type, ...
                                wavelength, medium, R0, d, N_inc, N_sca)
%simulate_shell Simulates the optical response of a (void) shell cluster
%
% High-level function to wrap a full simulation of core-shell system
% (to be called by a script with varying parameters)
%
% PARAMETERS:
% - rho_dye: surface density of shell dipoles
% - tensor_type: type of tensor [1 'isotropic', 2 'radial', 3 'random-flat', 4 'flat isotropic']
% - position_type: type of coverage [1 'fibonacci', 2 'random', >2 'hc']
% - wavelength: vector of wavelengths
% - medium: refractive index of surrounding medium
% - R0: sphere radius
% - d: dipole-sphere distance (for consistency with core-shell configuration)
% - N_inc: [optional, default: 36] number of incident angles
% - N_sca: [optional, default: 36] number of scattering angles
%
% RETURNS: matrix of cross-section spectra
%
% DEPENDS: alpha_bare, cluster_shell, spectrum_oa
%
% FAMILY: user_level, simulation
%

switch nargin
    case 7
        N_inc=36;
        N_sca=36;
    case 8
        N_sca=36;
end

N_dip = dye_coverage(rho_dye, R0);
if(mod(N_dip, 2) == 0)
    N_dip = N_dip+1;
end
N_dip

if(N_dip > 5e4)
    warning('probably dont want to do this (over 50k dipoles)')
end

switch position_type
    case 1 % fibonacci
        position = 'fibonacci';
    case 2 % random
        position = 'random';
    otherwise
        position = 'hc';
end

switch tensor_type
    case 1 % isotropic
        a = 1; b=1; c=1;
        orientation = 'radial'; %irrelevant
    case 2 % radial
        a = 0; b=0; c=1;
        orientation = 'radial'; %radial
    case 3 % random flat
        a = 0; b=0; c=1;
        orientation = 'flat'; %radial
    case 4 % flat-isotropic
        a = 1; b=1; c=0;
        orientation = 'radial'; %radial
    otherwise
        warning('Unexpected tensor type.')
end

alpha_inf=9.6e-39;
alpha_k = [5.76e-38];
lambda_k = [526];
mu_k = [10000];
exclusion = 0.5;
alphabar = alpha_bare(wavelength, alpha_inf, alpha_k, lambda_k, mu_k);
cl = cluster_shell(N_dip, R0, d, a, b, c, orientation, position, exclusion);

material.wavelength = wavelength(:);
material.alpha = alphabar;
% cater for dispersive media
if isnumeric(medium) && length(medium) == 1
    material.medium = material.wavelength*0 + medium;
else % all(ishandle(medium)) || all(ischar(medium))
    material.medium = feval(medium, wavelength);
end
xsec = spectrum_oa(cl, material, 'gl', N_inc, N_sca);
end
