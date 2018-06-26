function Rot = rotation_sphere(alpha, beta)
% rotation matrix (active)
% Note: equivalent to do
% Rot = rotation_euler_active(Alpha(j),Beta(j),Gamma(j));
% Rot = rotation_sphere(Alpha(j),Beta(j));
%
% PARAMETERS:
% - alpha: azimuth 
% - beta: polar angle
%
% RETURNS:
%
% DEPENDS: none
%
% FAMILY: low_level, utility, rotation
%
cosalpha = cos(alpha);  cosbeta = cos(beta);
sinalpha = sin(alpha);  sinbeta = sin(beta);

Rz = [cosalpha -sinalpha  0;
    sinalpha cosalpha  0;
    0      0       1]	;

Ry = [cosbeta 0  sinbeta;
    0       1        0;
    -sinbeta 0  cosbeta]	;

Rot = Rz*Ry;

end
