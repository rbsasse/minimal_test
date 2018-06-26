function Rot = rotation_euler_active(phi, theta, psi)
% Euler rotation matrix
% following the "zyz-convention"
% ie rotate(phi) around z, rotate(theta) around y', rotate(psi) around z''
% Note that this is an active rotation, ie to rotate colvec V, do R V,
% or for a row vector W=V^t, W R^t
%
% PARAMETERS:
% - phi: azimuth 
% - theta: polar angle
% - psi: polarisation angle
% 
% RETURNS:
%
% DEPENDS: none
%
% FAMILY: low_level, utility, rotation
%

cosphi = cos(phi); cospsi = cos(psi); costheta = cos(theta);
sinphi = sin(phi); sinpsi = sin(psi); sintheta = sin(theta);

Rot=zeros(3,3);

% note the indices below have been transposed from original code
Rot(1,1) = cosphi*costheta*cospsi - sinphi*sinpsi;
Rot(2,1) = sinphi*costheta*cospsi + cosphi*sinpsi;
Rot(3,1) = -sintheta*cospsi;

Rot(1,2) = -cosphi*costheta*sinpsi - sinphi*cospsi;
Rot(2,2) = -sinphi*costheta*sinpsi + cosphi*cospsi;
Rot(3,2) = sintheta*sinpsi;

Rot(1,3) = cosphi*sintheta;
Rot(2,3) = sinphi*sintheta;
Rot(3,3) = costheta;


end
