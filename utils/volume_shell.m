function V = volume_shell(radius, thickness)
%volume_shell Convert a surface coverage to a number of dyes
%
% PARAMETERS:
% - radius: surface density in nm^-2
% - thickness: radius of the shell
%
% RETURNS: number of dyes
%
% DEPENDS: none
%
% FAMILY: user_level, utility
%

  V = 4*pi/3*((radius + thickness)^3 - radius^3);

end
