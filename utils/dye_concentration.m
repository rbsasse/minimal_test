function [rho, N_dip] = dye_concentration(mD, mC, R)
%dye_concentration Convert a dye concentration to a surface density
%
% PARAMETERS:
% - mD: dye concentration in mol/L
% - mC: colloid concentration in mol/L
% - R: radius of the shell
%
% RETURNS: number of dyes
% - rho: surface density in nm^-2
%
% DEPENDS: none
%
% FAMILY: user_level, utility
%

if(nargin < 2)
    mC = 7e-12; % 7pM solution
    R = 30;
end

if(nargin < 3)
    R = 30;
end

Na = 6.022140857e23;
d_number = Na*mD; % number of dyes / L
c_number = Na*mC; % number of colloids / L
c_area = 4*pi*R^2;
total_area = c_number * c_area; % area in nm / L

rho = d_number / total_area;

N_dip = ceil(c_area * rho);

end

