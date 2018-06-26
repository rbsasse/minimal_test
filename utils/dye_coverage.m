function N_dip = dye_coverage( rho, R )
%DYE_COVERAGE Convert a surface coverage to a number of dyes
%
% PARAMETERS:
% - rho: surface density in nm^-2
% - R: radius of the shell
%
% RETURNS: number of dyes
%
% DEPENDS: none
%
% FAMILY: user_level, utility
%

area = 4*pi*R^2;
N_dip = ceil(area * rho);

end

