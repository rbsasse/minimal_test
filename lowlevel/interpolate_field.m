function Eint = interpolate_field(Egrid, theta)
%INTERPOLATE_FIELD interpolate the field components returned by Mie at arbitrary thetas
%
% PARAMETERS:
% - Egrid Structure with fields theta, Ecr, Ect, Esf, [Em0r, Em0t] (dipole excitation)
% - theta vector of angles to interpolate at
%
% RETURNS: Structure with fields Ecr, Ect, Esf, [Em0r, Em0t] (dipole excitation)
%
% DEPENDS: none
%
% FAMILY: low_level, utility
%

thetagrid = real(Egrid.theta);
Eint.Ecr  = interp1(thetagrid, real(Egrid.Ecr), theta) + ...
    1i*interp1(thetagrid, imag(Egrid.Ecr), theta);
Eint.Ect  = interp1(thetagrid, real(Egrid.Ect), theta) + ...
    1i*interp1(thetagrid, imag(Egrid.Ect), theta);
Eint.Esf  = interp1(thetagrid, real(Egrid.Esf), theta) + ...
    1i*interp1(thetagrid, imag(Egrid.Esf), theta);

% these two fields don't exist for PWE
if isfield(Egrid,'Em0r')
    Eint.Em0r = interp1(thetagrid, real(Egrid.Em0r), theta) + ...
        1i*interp1(thetagrid, imag(Egrid.Em0r), theta);
end
if isfield(Egrid,'Em0t')
    Eint.Em0t = interp1(thetagrid, real(Egrid.Em0t), theta) + ...
        1i*interp1(thetagrid, imag(Egrid.Em0t), theta);
end

end
