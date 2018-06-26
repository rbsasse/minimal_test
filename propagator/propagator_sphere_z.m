function hSc = propagator_sphere_z(theta, phi, Egrid)
%PROPAGATOR_SPHERE_Z Sphere's Green tensor in the rotated frame (source along z)

% each column corresponds to the field for a unit dipole
% located along the z-axis
% and with three orthogonal orientations (x, y, z)
% Sgrid a 3nThetax3 matrix (nTheta 3x3 blocks) to interpolate
%
% PARAMETERS:
% - theta
% - phi
% - Egrid
%
% RETURNS: 3Nrx3 matrix, each 3x3 block being the Green's tensor in cartesian coordinates
%
% DEPENDS: interpolate_field
%
% FAMILY: low_level, propagator, sphere
%

%% convenience
cosphi = cos(phi);
sinphi = sin(phi);
costheta = cos(theta);
sintheta = sin(theta);

%% interpolated field components
Edint = interpolate_field(Egrid, theta);

%% re-including phi dependence
%
% field due to dipole along x
pxEr = Edint.Ecr.*cosphi;    % radial component
pxEt = Edint.Ect.*cosphi;    % theta  component
pxEf = Edint.Esf.*sinphi;    % phi component

% field due to dipole along y
% phi -> phi - pi/2 
% hence cos(phi) -> sin(phi), sin(phi) -> -cos(phi)
pyEr = Edint.Ecr.*sinphi;    % radial component
pyEt = Edint.Ect.*sinphi;           % theta  component
pyEf = -1*Edint.Esf.*cosphi;           % phi component

% field due to dipole along z
pzEr = Edint.Em0r;    % radial component
pzEt = Edint.Em0t;    % theta  component
% pzEf = 0*theta; % no phi component

%% transform fields to cartesian
% note: the idea is to rotate, but can't use rotation matrix
% since each row requires its own rotation

% NOTE: should use Ec(px) = spherical_to_cartesian(theta, phi, Es(px))
pxEx = pxEr.*sintheta.*cosphi + pxEt.*costheta.*cosphi - pxEf.*sinphi;
pxEy = pxEr.*sintheta.*sinphi + pxEt.*costheta.*sinphi + pxEf.*cosphi;
pxEz = pxEr.*costheta         - pxEt.*sintheta                       ;

pyEx = pyEr.*sintheta.*cosphi + pyEt.*costheta.*cosphi - pyEf.*sinphi;
pyEy = pyEr.*sintheta.*sinphi + pyEt.*costheta.*sinphi + pyEf.*cosphi;
pyEz = pyEr.*costheta         - pyEt.*sintheta                       ;

pzEx = pzEr.*sintheta.*cosphi + pzEt.*costheta.*cosphi               ;
pzEy = pzEr.*sintheta.*sinphi + pzEt.*costheta.*sinphi               ;
pzEz = pzEr.*costheta         - pzEt.*sintheta                       ;

% bind as three rows
%Ex1 Ex2 Ex3 Ex4 Ex5 ...
%Ey1 Ey2 Ey3 Ey4 Ey5 ...
%Ez1 Ez2 Ez3 Ez4 Ez5 ...

pxEc = [pxEx; pxEy; pxEz]; % First column of S' for each dipole
pyEc = [pyEx; pyEy; pyEz];
pzEc = [pzEx; pzEy; pzEz];

%% returning a 3Nx3 matrix
% each 3x3 block is the Green's tensor in cartesian coordinates
% i.e first  column is Ex, Ey, Ez for dipole along x
%     second column is Ex, Ey, Ez for dipole along y
%     third  column is Ex, Ey, Ez for dipole along z
% the prefactor Ep0 is according to the Mie definition
hSc = Egrid.Ep0*[pxEc(:), pyEc(:), pzEc(:)]; % [3N x 3]
% TODO check that the blocks aren't somehow transposed
end
