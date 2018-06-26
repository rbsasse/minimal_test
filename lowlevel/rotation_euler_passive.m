function Rot = rotation_euler_passive(phi, theta, psi)
% Euler rotation matrix
% following the "zyz-convention"
% ie rotate(phi) around z, rotate(theta) around y', rotate(psi) around z''
% Note that this is a passive rotation (ie of the frame),
% ie to rotate colvec V, do R^t V,
% or for a row vector W=V^t, W R
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

Rot(1,1) = cosphi*costheta*cospsi - sinphi*sinpsi;
Rot(1,2) = sinphi*costheta*cospsi + cosphi*sinpsi;
Rot(1,3) = -sintheta*cospsi;

Rot(2,1) = -cosphi*costheta*sinpsi - sinphi*cospsi;
Rot(2,2) = -sinphi*costheta*sinpsi + cosphi*cospsi;
Rot(2,3) = sintheta*sinpsi;

Rot(3,1) = cosphi*sintheta;
Rot(3,2) = sinphi*sintheta;
Rot(3,3) = costheta;


end
