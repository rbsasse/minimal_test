function  xsec = simulate_mie(rho_dye, R0, d, spacer, wavelength, medium, ...
      N_nmax, shell, mloc)
  %simulate_mie Simulates the optical response of homogeneous Mie shell
  %
  % High-level function to wrap a full simulation of homogeneous shell system
  % (to be called by a script with varying parameters)
  %
  % PARAMETERS:
  % - rho_dye: surface density of shell dipoles
  % - wavelength: vector of wavelengths
  % - medium: refractive index of surrounding medium
  % - R0: sphere radius
  % - d: thickness of dye layer
  % - spacer: thickness of spacer layer
  % - N_nmax: [optional, default: 20] number of incident angles
  % - shell: logical, void or silver core
  % - mloc: logical, normalise by mloc
  %
  % RETURNS: matrix of cross-section spectra
  %
  % DEPENDS: epsilon_Ag, alpha_bare, volume_shell, epsilon_dye, mini_FF
  %
  % FAMILY: user_level, simulation
  %

switch nargin
    case 8
        spacer=1;
    case 7
        spacer=1;
        mloc=false;
    case 6
        spacer=1;
        shell = false;
        mloc=false;
    case 5
        spacer=1;
        N_nmax = 20;
        shell = false;
        mloc=false;
end

wavelength=wavelength(:);
% cater for dispersive media
if isnumeric(medium) && length(medium) == 1
    medium = wavelength*0 + medium;
else % all(ishandle(medium)) || all(ischar(medium))
    medium = feval(medium, wavelength);
end

Ca={R0, R0+spacer, R0+spacer+d};

% This corresponds to alpha_{zz} for R6G assuming a single Lorentzian
alpha_inf=9.6e-39;
alpha_k = 5.76e-38;
lambda_k = 526;
mu_k = 10000;
eps0 = 8.854e-12;
nm3 = 1e27;
prefact = nm3/(4*pi*eps0); % scaling in cd code, to invert here
alpha = 1/prefact * alpha_bare(wavelength, alpha_inf, alpha_k, lambda_k, mu_k);

shell_surface = 4*pi*(R0+spacer+d)^2;
N_dye = rho_dye * shell_surface

Vshell = volume_shell(R0, d);
c_dye = N_dye / Vshell; % dye/nm^3 in shell
eps_dye = epsilon_dye(alpha, c_dye*nm3, medium);
eps_medium = medium.^2;
if(shell)
    eps_core = eps_medium;
else
    eps_core = epsilon_Ag(wavelength);
end
Cepsilon={eps_core, eps_medium, eps_dye, eps_medium};
Cepsilon2={eps_core, eps_medium};
coated = mini_FF(wavelength,Cepsilon,Ca,N_nmax);
bare = mini_FF(wavelength,Cepsilon2,{R0},N_nmax);
%
% nNbTheta = 300;
% mid = R0+d;
% stE = mini_NF(coated,mid,nNbTheta);

xsec = ([coated.Qext coated.Qabs coated.Qsca] * pi * (R0+d+spacer)^2 - ...
    [bare.Qext bare.Qabs bare.Qsca] * pi * (R0)^2 ) / N_dye ;
if(mloc)
    xsec = bsxfun(@times, xsec, 1./coated.MLocAve);
end
end
