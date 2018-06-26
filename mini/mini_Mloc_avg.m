function MLocAve = mini_Mloc_avg(x,stGD)
% Calculates the surface-averaged local field intensity EF for PWE
% from the Mie susceptibilities using Eq. (H.79).
%
% Parameters:
% - x:      column vector [L x 1]
%           wavelength-dependent x=kM*a (Eq. H.45)
% - stGD:  structure with two fields, Gamma and Delta
%           each a matrix [L x nNmax]
%           with suceptibilities Gamma_n and Delta_n
%           can be obtained from function GenSuscepGDAB
%
% Returns:
% - MLocAve: [L x 1] wavelength-dependent average LFIEF
%
% This file is part of the SPlaC v1.0 package (copyright 2008)
% Check the README file for further information

nNmax=size(stGD.Gamma,2);

% get psi_n(x), xi_n(x) and derivatives
stRBx=mini_RB(nNmax,x);

n=transpose(1:nNmax); % [nNmax x 1]

cc1 = 2.*n+1; % [nNmax x 1]
% Get sums and check loss of precision during sums
tmp1=stRBx.psi + stGD.Gamma .* stRBx.xi;
tmp2=stRBx.Dpsi + stGD.Delta .* stRBx.Dxi;
summat1=abs(tmp1).^2 + abs(tmp2).^2; % [L x nNmax]

clear tmp2;
cc2 = cc1.*n.*(n+1); % [nNmax x 1]
tmp1=stRBx.psi + stGD.Delta .* stRBx.xi;
summat2=abs(tmp1).^2; % [L x nNmax]

% From Eq. H.79
MLocAve= 1./(2*(x).^2) .* (summat1 * cc1) ...
    + 1./(2*(x).^4) .* (summat2 * cc2) ; % [L x 1]
% sums over positive numbers: no loss of precision issue

end
