function Egrid = mie_PWE_z(d, wavelength, epsilon, medium, radius, ...
    N_nmax, N_theta,stPiTauMie)
%mie_PWE_Z Evalutes the sphere's scattered field for PWE illumination along z
%
% The scattered field is evaluated for two orthogonal polarisations
% and at multiple angles theta in the local frame. These values will be
% interpolated for actual dipole positions in the lab frame in incident_field_sphere
%
% PARAMETERS:
% - d distance from core to dipoles
% - wavelength
% - epsilon dielectric function of core
% - medium refractive index of surrounding medium
% - radius radius of the core sphere
% - N_nmax number of multipoles considered in the solution
% - N_theta number of theta points to evalute the fields at
% - stPiTauMie [optional] precalculated angular functions pi and tau
%
% RETURNS: structure with grid of field components under PWE along z
%
% DEPENDS: mie_PWE_Etheta
%
% FAMILY: low_level, mie
%

global noCheckSum;
noCheckSum=true;


theta=transpose(linspace(0,pi,N_theta)); % row [1 x T]
k0 = 2*pi/wavelength;
kM = k0*medium;
s = sqrt(epsilon) / medium;

if nargin<10
    stPiTauMie = mie_vshPinmTaunm01(N_nmax,theta);
end


stParam.a=radius;
stParam.d=d;
stParam.s=s;
stParam.x=kM*radius;

r = radius + d;

stE = mie_PWE_Etheta(stParam, N_nmax, theta, r, stPiTauMie);

Egrid.theta = theta;
Egrid.Ecr = stE.Ecr;
Egrid.Ect = stE.Ect;
Egrid.Esf = stE.Esf;

end
