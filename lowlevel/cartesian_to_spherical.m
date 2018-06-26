function [ Vs ] = cartesian_to_spherical(theta, phi, Vc)
%% cartesian_to_spherical
% Convert vector components from spherical to
% cartesian coordinates, for N vectors at a time
% 
%	cartesian_to_spherical(theta, phi, Vc) converts Vc to spherical
%
% Input:
% - theta: [1 x N_dip] vector of polar angles
% - phi:   [1 x N_dip] vector of azimutal angles
% - Vc: complex [3 x N_dip] matrix of cartesian vector components (column-wise)
%
% Output: 
% - Vs: complex [3 x N_dip] matrix of spherical vector components (column-wise)
%
% Note: the idea is to rotate, but can't use rotation matrix
% since each row requires its own rotation
%
% Dependency: 
% none
phi = phi(:)'; theta = theta(:)'; % ensure row vectors

cosphi = cos(phi);
sinphi = sin(phi);
costheta = cos(theta);
sintheta = sin(theta);

% V_r
Vs(1,:) =  Vc(1,:).*sintheta.*cosphi + Vc(2,:).*sintheta.*sinphi + Vc(3,:).*costheta;
% V_theta
Vs(2,:) =  Vc(1,:).*costheta.*cosphi + Vc(2,:).*costheta.*sinphi - Vc(3,:).*sintheta;
% V_phi
Vs(3,:) = -Vc(1,:).*sinphi           + Vc(2,:).*cosphi                              ;

end

