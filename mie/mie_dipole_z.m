function Egrid = mie_dipole_z(d, wavelength, epsilon, medium, radius, ...
    N_esa, N_mie, N_theta, stPiTauMie,stPiTauESA)
%mie_DIPOLE_Z Wrapper for the Mie solution of dipole excitation along z axis
%
% Evaluates the field factors for the sphere's Green tensor
% at many angles theta in the local frame where the source is on the z axis.
% This function is merely a wrapper for mie_DipEthetaFromMieAndESA with more convenient
% parameters for the problem at hand.
% FUTURE: Will also return the far-field at grid of angles, so that
% one may compute FF cross-sections of the sphere due to excitation by dipoles
%
% PARAMETERS:
% - stParam structure of parameters, incl. a, d, s, x
% - N_esa maximum multipole order considered in total
% - N_mie maximum multipole order considered exactly
% - theta polar angles where the field is to be evaluated
% - r radius where the field is to be evaluated
% - stPiTau [optional] precomputed pi and tau angular functions
%
% RETURNS: a structure with the field components at grid of positions
%
% DEPENDS: mie_vshPinmTaunm01, mie_DipEthetaFromMieAndESA
%
% FAMILY: low_level, mie
%

global noCheckSum;
noCheckSum=true;


theta=transpose(linspace(0,pi,N_theta)); % row [1 x T]
k0 = 2*pi/wavelength;
kM = k0*medium;
s = sqrt(epsilon) / medium;

if nargin < 9
    stPiTauMie = mie_vshPinmTaunm01(N_mie,theta);
    stPiTauESA = mie_vshPinmTaunm01(N_esa,theta);
end


stParam.a=radius;
stParam.d=d;
stParam.s=s;
stParam.x=kM*radius;

r = radius + d; % we happen to need the field on the same shell as the source

stE = mie_DipEthetaFromMieAndESA(stParam,N_mie,N_esa,theta,r,stPiTauMie,stPiTauESA);

% Eq. H83
% prefactor for dipole VSWF expansion in medium
% note: |p| introduced elsewhere
% note: \bar p = p / kappa -> 4pi epsilon_0 epsilon_m factor
% note: sqrt(6pi) swallowed by our definition of a's b's w.r.t. [H.84]
Egrid.Ep0 = 4*pi*1i*kM.^3 ;
Egrid.theta = theta;
% field components - note 'c' stands for cos(phi), 's' for sin(phi) factors
% to be re-introduced when calculating fields at arbitrary positions
Egrid.Em0r = stE.Em0r; % radial component for perp dipole (pz)
Egrid.Em0t = stE.Em0t; % theta component for perp dipole (pz)
Egrid.Ecr = stE.Ecr; % radial component for parallel dipole (px)
Egrid.Ect = stE.Ect; % theta component for parallel dipole (px)
Egrid.Esf = stE.Esf; % phi component for parallel dipole (px)

% TODO return also FF components

end
