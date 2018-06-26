function stE = mini_NF(stM,r0,nNbTheta)

% Calculates several properties of the field on a sphere surface at r=r0 for PWE
% For PWE, we have: |m|=1; a,c even in m; b,d, odd in m.
% The 4 Mie coefficients in stAabcdn1 must be the relevant ones
% at r=r0.
% Use r0=Inf to obtain far-field properties.
% The fields Ecr, Ect, Esf given in the results are discussed in the
% supplementary information.
%
% Parameters:
% - lambda: column vector [L x 1]
%           wavelengths in nm
% - epsilon: scalar or column vector [L x 1]
%           epsilon of dielectric where field is evaluated
% - stAbcdn1: structure with 4 fields, an1, bn1, cn1, dn1,
%             each a matrix [L x nNmax]
%             containing the Mie coefficients a_{n,1}, b_{n,1},
%             c_{n,1} and d_{n,1} used for the VSH expansions at r=r0.
% - r0:     scalar [1 x 1] non-zero
%           distance from origin (in nm)
%           r0 can be "Inf" to obtain far-field radiation profile
% - nNbTheta: integer scalar
%             number of theta points used for computations
% - sRegion: string
%            'outside', 'inside', 'scattering, or 'all' (default):
%            * if 'outside': the points should be in the outside region (outside the largest
%              sphere). The regular part of the field (coeffs a and b) is then known and not
%              computed from the plane wave incident field (faster).
%            * if 'scattering': the regular part of the field (coeffs a and
%              b) is set to zero (no a and b needed)
%            * if 'inside': the points should be in the inside region (containing the origin), then
%              the irregular part of the field is zero (no coeffs c and d needed).
%            * otherwise: all of a,b,c,d are used.
%
% Returns: structure stEsurf with fields
% - stEsurf.theta: [1 x T] row vector with theta's
% - stEsurf.r0: [1 x 1] r0 (in nm)
% - stEsurf.lambda: [L x 1] lambda (in nm)
% - stEsurf.Ecr: [L x T] wavelength-and-theta-dependent Ecr
% - stEsurf.Ect: [L x T] wavelength-and-theta-dependent Ect
% - stEsurf.Esf: [L x T] wavelength-and-theta-dependent Esf
% - stEsurf.MLocParaAve: [L x 1] wavelength-dependent average MLocPara
% - stEsurf.MLocPerpAve: [L x 1] wavelength-dependent average MLocPerp
% - stEsurf.MLocAve: [L x 1] wavelength-dependent average MLoc
% - stEsurf.F0E4Ave: [L x 1] wavelength-dependent average SERS EF F^0_{E4}
% - stEsurf.F0E4PerpAve: [L x 1] wavelength-dependent average SERS EF
%                        F^0_{E4} for perpendicular component
% - stEsurf.F0E4ParaAve: [L x 1] wavelength-dependent average SERS EF
%                        F^0_{E4} for parallel component
%
% This file is part of the SPlaC v1.0 package (copyright 2008)
% Check the README file for further information


if (nNbTheta<91)
    disp 'Warning in PweEsurf: nNbTheta must be large for averages to be meaningful'
end

theta=linspace(0,pi,nNbTheta); % row [1 x T]
stPinTaun = mini_pitau(stM.nNmax,theta(:));

% stM corresponds to a spherical multilayer
% find layer where r0 is located
nK=stM.nK;
kk=nK;
while ((kk>0) && (r0<0.99999*stM.Ca{kk})) % note that if r0=Ca{kk}, properties are computed just outside
    kk=kk-1;
end
% r0 has been located in region kk
% extract corresponding Mie coefficients
stAbcdn1.an1=stM.Calphakn1{kk+1};
stAbcdn1.bn1=stM.Cbetakn1{kk+1};
stAbcdn1.cn1=stM.Cgammakn1{kk+1};
stAbcdn1.dn1=stM.Cdeltakn1{kk+1};

stE_reg=mini_E_PWE(stM.lambda,stM.Cepsilon{kk+1},stAbcdn1.an1,stAbcdn1.bn1,r0,theta(:),'j',stPinTaun);
stE_irr=mini_E_PWE(stM.lambda,stM.Cepsilon{kk+1},stAbcdn1.cn1,stAbcdn1.dn1,r0,theta(:),'h1',stPinTaun);

% disp 'PweEsurf: Compiling results and averages...'
% Ecr, Ect, and Esf are [L x T]
stE.Ecr= stE_reg.Ecr + stE_irr.Ecr;
stE.Ect= stE_reg.Ect + stE_irr.Ect;
stE.Esf= stE_reg.Esf + stE_irr.Esf;

stE.theta=theta;
stE.r0=r0;
stE.lambda=stM.lambda;

clear stE_reg stE_irr;

% computes surface-averages

Ecr2theta=abs(stE.Ecr).^2;
Ect2theta=abs(stE.Ect).^2;
Esf2theta=abs(stE.Esf).^2;

dtheta=pi/(nNbTheta-1);
sintcolnorm=(dtheta/2)*transpose(sin(theta)); % [T x 1]

% average LFIEF
stE.MLocPerpAve=1/2* Ecr2theta * sintcolnorm;
stE.MLocParaAve=1/2* (Ect2theta +Esf2theta) * sintcolnorm;
stE.MLocAve=stE.MLocPerpAve+stE.MLocParaAve;


end
