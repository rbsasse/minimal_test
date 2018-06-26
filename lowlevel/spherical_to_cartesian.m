function [ Vc ] = spherical_to_cartesian(theta, phi, Vs)
%% spherical_to_cartesian
% Convert vector components from spherical to
% cartesian coordinates, for multiple vectors at a time
% 
%	spherical_to_cartesian(theta, phi, Vs) converts Vs to cartesian
%
% Input:
% - theta: [1 x N] vector of polar angles
% - phi:   [1 x N] vector of azimutal angles
% - Vs: complex [3 x Nr] matrix of spherical vector components (column-wise)
%
% Output: 
% - Vc: complex [3 x Nr] matrix of cartesian vector components (column-wise)
%
% Note: the idea is to rotate, but can't use rotation matrix
% since each row requires its own rotation
%
% Dependency: 
% none

cosphi = cos(phi);
sinphi = sin(phi);
costheta = cos(theta);
sintheta = sin(theta);

% V_x
Vc(1,:) = Vs(1,:).*sintheta.*cosphi + Vs(2,:).*costheta.*cosphi - Vs(3,:).*sinphi;
% V_y
Vc(2,:) = Vs(1,:).*sintheta.*sinphi + Vs(2,:).*costheta.*sinphi + Vs(3,:).*cosphi;
% V_z
Vc(3,:) = Vs(1,:).*costheta         - Vs(2,:).*sintheta                          ;


end

