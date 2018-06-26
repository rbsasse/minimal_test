function [stE, stEn] = mie_DipEthetaFromMie(stParam,N,theta,r,stPiTau)
  %mie_DIPETHETAFROMMIE Dipolar field scattered by a sphere (Mie solution)
  %
  % Mie solution for the scattered field from a dipole along the z-axis in the presence of a sphere
  %
  % PARAMETERS:
  % - stParam structure of parameters, incl. a, d, s, x
  % - N maximum multipole order considered
  % - theta polar angles where the field is to be evaluated
  % - r radius where the field is to be evaluated
  % - stPiTau [optional] precomputed pi and tau angular functions
  %
  % RETURNS: a structure with the field components at grid of positions
  %
  % DEPENDS: mie_vshPinmTaunm01, RBpsi2 [private], RB [private], ZnAllIrr [private]
  %
  % FAMILY: low_level, mie
  %

  if nargin<5
    stPiTau = mie_vshPinmTaunm01(N,theta);
  end

  a=stParam.a;
  d=stParam.d;
  % lambda=stParam.lambda;
  % epsilonM=stParam.epsilonM;
  s=stParam.s;
  x=stParam.x;
  %kM=2*pi/lambda*sqrt(epsilonM);
  %x=kM*a;
  xp=x/a*(a+d);

  nvec=1:N; % [1 x N]
  z=s.*x;

  stRBz=RBpsi2(N,z);
  stRBx=RB(N,x);

  scale=1e150;
  invsc=1e-150;

  PP1 = (scale*stRBz.psi) .* stRBx.Dpsi;
  PP2 = (scale*stRBx.psi) .* stRBz.Dpsi;
  PP3 = stRBz.psi .* stRBx.Dxi;
  PP4 = stRBx.xi .* stRBz.Dpsi;

  % numerators
  NumGam = - PP1 + s * PP2;
  NumDel = PP2 - s * PP1;

  % denominators
  DenDelB = - PP4 + s * PP3;
  DenGamA = PP3 - s * PP4;

  stZxp=ZnAllIrr(N,xp);
  % Z0 is xi(xp)/xp, Z2 is (x xi)'(xp)/xp
  % Z1 is xi(xp)/xp^2

  % Parallel dipole
  coeff1=sqrt(3/8*(2*nvec+1))*invsc; % [1 x N]

  an1 = 1i * coeff1 .* stZxp.Z0;
  bn1 = - coeff1 .* stZxp.Z2;

  % Perpendicular dipole
  coeff2=2*coeff1.*sqrt(nvec.*(nvec+1)); % [1 x N]

  bn0 = coeff2 .* stZxp.Z1;

  % calculates susceptibilities from Eqs. H.46, H.47, H. 51, H.52
  % Gamman = NumGam ./ (DenGamA*scale);
  % Deltan = NumDel ./ (DenDelB*scale);

  % calculates scattered field coeffs
  cn1=(an1 .* NumGam) ./ DenGamA;
  dn1=(bn1 .* NumDel) ./ DenDelB;
  dn0=(bn0 .* NumDel) ./ DenDelB;

  stZnAll=ZnAllIrr(N,x*r/a); % fields are [1 x N]

  nnp1vec=nvec.*(nvec+1);
  kappan=sqrt((2*nvec+1)./(4*pi*nnp1vec));
  dwign1=bsxfun(@times,stPiTau.pin1,sin(theta));

  if nargout==1 % Fast version, no n-dependence returned

    tmpnr=nnp1vec.*kappan.*stZnAll.Z1;
    tmp=tmpnr.*dn0;
    stE.Em0r=stPiTau.Pn0 * tmp.'; % [T x 1]

    tmpnt=kappan.*stZnAll.Z2;
    tmp=tmpnt.*dn0;
    stE.Em0t=stPiTau.taun0 * tmp.'; % [T x 1]

    tmp=-2*tmpnr.*dn1;
    stE.Ecr=dwign1 * tmp.'; % [T x 1]

    tmp=-2*tmpnt.*dn1;
    EctN=stPiTau.taun1 * (tmp).'; % [T x 1]
    EsfN=stPiTau.pin1 * (-tmp).'; % [T x 1]

    tmp=-2i*kappan.*stZnAll.Z0.*cn1;
    EctM=stPiTau.pin1 * (tmp).'; % [T x 1]
    EsfM=stPiTau.taun1 * (-tmp).'; % [T x 1]

    stE.Ect=EctN+EctM;
    stE.Esf=EsfN+EsfM;

  else % Slow version with all n-terms in series

    tmpnr=nnp1vec.*kappan.*stZnAll.Z1;
    tmp=tmpnr.*dn0;
    stEn.Em0rn=bsxfun(@times,stPiTau.Pn0,tmp); % [T x N]
    stE.Em0r=stPiTau.Pn0 * tmp.'; % [T x 1]

    tmpnt=kappan.*stZnAll.Z2;
    tmp=tmpnt.*dn0;
    stEn.Em0tn=bsxfun(@times,stPiTau.taun0,tmp); % [T x N]
    stE.Em0t=stPiTau.taun0 * tmp.'; % [T x 1]

    tmp=-2*tmpnr.*dn1;
    stEn.Ecrn=bsxfun(@times,dwign1,tmp); % [T x N]
    stE.Ecr=dwign1 * tmp.'; % [T x 1]

    tmp=-2*tmpnt.*dn1;
    stEn.EctNn=bsxfun(@times,stPiTau.taun1,tmp); % [T x N]
    EctN=stPiTau.taun1 * (tmp).'; % [T x 1]
    stEn.EsfNn=bsxfun(@times,stPiTau.pin1,-tmp); % [T x N]
    EsfN=stPiTau.pin1 * (-tmp).'; % [T x 1]

    tmp=-2i*kappan.*stZnAll.Z0.*cn1;
    stEn.EctMn=bsxfun(@times,stPiTau.pin1,tmp); % [T x N]
    EctM=stPiTau.pin1 * (tmp).'; % [T x 1]
    stEn.EsfMn=bsxfun(@times,stPiTau.taun1,-tmp); % [T x N]
    EsfM=stPiTau.taun1 * (-tmp).'; % [T x 1]

    stEn.Ectn = stEn.EctNn + stEn.EctMn;
    stEn.Esfn = stEn.EsfNn + stEn.EsfMn;
    stE.Ect=EctN+EctM;
    stE.Esf=EsfN+EsfM;

  end
end

function stRBpsi2=RBpsi2(nNmax, rho)
  % rho scalar here

  n=1:nNmax;
  nm1=0:nNmax;
  nu=nm1+0.5;

  f = besselj(nu,rho);
  % f is matrix [1 x nNmax+1] of cylindrical Bessel
  % J_{n+0.5}(rho), n=0..nNmax

  sq=sqrt((pi/2)*rho); % [1 x 1]
  f=sq*f; % [1 x nNmax+1]
  % f is now row of spherical Ricati-Bessel
  % psi_n(rho), n=0..nNmax or equivalently psi_{n-1}(rho), n=1..nNmax+1

  stRBpsi2.psi=f(1,n+1);

  % Computes: psi_n'=psi_{n-1} - n psi_n/rho
  % and check for loss of precision in sum
  stRBpsi2.Dpsi=f(1,n) - ( (1/rho)*n ).*f(1,n+1);
end

function stRB=RB(nNmax, rho)
  % rho scalar here

  n=1:nNmax;
  nm1=0:nNmax;
  nu=nm1+0.5;

  fj=besselj(nu,rho);
  f=besselh(nu,rho);

  sq=sqrt((pi/2)*rho); % [1 x 1]
  fj=sq*fj; % [R x nNmax+1]
  % fj is now matrix of spherical Ricati-Bessel
  % psi_n(rho), n=0..nNmax or equivalently psi_{n-1}(rho), n=1..nNmax+1
  f=sq*f; % [R x nNmax+1]
  % f is now matrix of Ricatti-Bessel
  % xi_n(rho), n=0..nNmax or equivalently xi_{n-1}(rho), n=1..nNmax+1

  % Computes: psi_n'=psi_{n-1} - n psi_n/rho
  % and check for loss of precision in sum
  stRB.Dpsi=fj(1,n) - ( (1/rho)*n ).*fj(1,n+1);

  % Computes: xi_n'=xi_{n-1} - n xi_n/rho
  % and check for loss of precision in sum
  stRB.Dxi=f(1,n) - ( (1/rho)*n ).*f(1,n+1);

  stRB.psi=fj(1,n+1);
  stRB.xi=f(1,n+1);

end

function stZnAll=ZnAllIrr(nNmax, rho)
  % rho scalar here

  n=1:nNmax;
  nm1=0:nNmax;
  nu=nm1+0.5;

  f=besselh(nu, rho);

  % f is matrix [1 x nNmax+1] of cylindrical Bessel
  % Z_{n+0.5}(rho), n=0..nNmax

  sq=sqrt((pi/2)/rho); % [1 x 1]
  f=sq*f; % [R x nNmax+1]
  % f is now matrix of spherical Bessel
  % z_n(rho), n=0..nNmax or equivalently z_{n-1}(rho), n=1..nNmax+1

  stZnAll.Z0=f(1,2:(nNmax+1));
  stZnAll.Z1=stZnAll.Z0/rho;

  % Computes: Z2_n=z_{n-1} - n Z1_n
  stZnAll.Z2=f(1,n) - n.*stZnAll.Z1;

end
