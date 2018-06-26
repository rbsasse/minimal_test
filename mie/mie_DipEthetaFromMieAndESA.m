function stE = mie_DipEthetaFromMieAndESA(stParam, N_mie, N_esa, theta, r,...
    stPiTauMie,stPiTauESA)
%mie_DIPETHETAFROMMIEANDESA Dipolar field scattered by a sphere (full solution)
%
% Combined Mie and ESA solutions for the scattered field from a dipole along the z-axis in the presence of a sphere.
% The Mie solution is accurate for small multipole orders, typically up to about 30-50;
% from there the ESA offers a more robust approach. We therefore add the two contributions,
% subtracting the lower-order terms from the ESA solution that have been computed rigorously with Mie.
%
% PARAMETERS:
% - stParam structure of parameters, incl. a, d, s, x
% - NMie maximum multipole order considered in the Mie solution
% - NESA maximum multipole order considered in the ESA solution
% - theta polar angles where the field is to be evaluated
% - r radius where the field is to be evaluated
% - stPiTauMie [optional] precomputed pi and tau angular functions
% - stPiTauESA [optional] precomputed pi and tau angular functions
%
% RETURNS: a structure with the field components at grid of positions
%
% DEPENDS: mie_shPinmTaunm01, mie_DipEthetaFromMie, mie_DipEthetaFromESA
%
% FAMILY: low_level, mie
%

if nargin<6
    stPiTauMie = mie_vshPinmTaunm01(N_mie,theta);
    stPiTauESA = mie_vshPinmTaunm01(N_esa,theta);
end

stM = mie_DipEthetaFromMie(stParam,N_mie,theta,r,stPiTauMie);
if N_esa == 0 % testing Mie only
    stE = stM;
else
    stEshort = mie_DipEthetaFromESA(stParam,N_mie,theta,r,stPiTauMie);
    
    stEfull = mie_DipEthetaFromESA(stParam,N_esa,theta,r,stPiTauESA);
    
    stE.Em0r = stM.Em0r + stEfull.Em0r - stEshort.Em0r;
    stE.Em0t = stM.Em0t + stEfull.Em0t - stEshort.Em0t;
    stE.Ecr = stM.Ecr + stEfull.Ecr - stEshort.Ecr;
    stE.Ect = stM.Ect + stEfull.Ect - stEshort.Ect;
    stE.Esf = stM.Esf + stEfull.Esf - stEshort.Esf;
end
end
