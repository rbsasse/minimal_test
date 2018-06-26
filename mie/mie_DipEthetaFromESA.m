function [stE, stEn] = mie_DipEthetaFromESA(stParam,N,theta,r,stPiTau)
  %mie_DipEthetaFromESA Dipolar field scattered by a sphere in the electrostatic approximation
  %
  % ESA part of the scattered field from a dipole along the z-axis in the presence of a sphere
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
  % DEPENDS: mie_vshPinmTaunm01
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
  R=a+d;

  nvec=1:N; % [1 x N]

  E0ovEP=-1i*sqrt(3/(8*pi))/x^3;

  geomf=(a^2/r/R).^(nvec+2);
  betan=(s^2-1)./(s^2+1+1./nvec);
  bgn=E0ovEP*betan.*geomf; % [1 x N]

  sbgn=sqrt(nvec.*(nvec+1)).*bgn; % [1 x N]

  if nargout==1 % Fast version, only sums

    stE.Em0r=stPiTau.Pn0 *((nvec+1).^2.*bgn).'; % [T x 1]
    stE.Em0t= stPiTau.taun0 * (-(nvec+1).*bgn).';

    stE.Ecr = sin(theta).* (stPiTau.pin1 * (-(nvec+1).*sbgn).' );
    stE.Ect = stPiTau.taun1 * sbgn.';
    stE.Esf = stPiTau.pin1 * (-sbgn).';

  else

    stEn.Em0rn = bsxfun(@times,stPiTau.Pn0,(nvec+1).^2.*bgn); % [T x N]
    stEn.Em0tn = bsxfun(@times,stPiTau.taun0,-(nvec+1).*bgn);

    stEn.Ecrn = ( sin(theta) * (-(nvec+1).*sbgn) ) .* stPiTau.pin1;
    stEn.Ectn = bsxfun(@times,stPiTau.taun1,sbgn);
    stEn.Esfn = bsxfun(@times,stPiTau.pin1,-sbgn);

    stE.Em0r=sum(stEn.Em0rn,2);
    stE.Em0t=sum(stEn.Em0tn,2);
    stE.Ecr=sum(stEn.Ecrn,2);
    stE.Ect=sum(stEn.Ectn,2);
    stE.Esf=sum(stEn.Esfn,2);
  end

end
