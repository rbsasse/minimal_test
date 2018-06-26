function stMieab = mie_ab_geom_PWE(theta, phi, N_nmax)
%% mie_ab_geom_PWE
% Calculates geometric part of Mie coefficients a and b for PWE along
% specific directions
%
%	mie_ab_geom_PWE(theta, phi, N_nmax) computes angular amn, bmn
%
% These coefficients need to be combined with the incident field
% spherical components to produce the actual a and b of the multipole expansion.
%
% Input:
% - theta: numeric vector [N_inc x 1] PW incidence angle
% - phi: numeric vector [N_inc x 1] PW incidence angle
% - N_nmax: scalar integer [1 x 1]
%
% Output: structure with fields:
% - a_theta: matrix [N_inc x nP]
% - a_phi: matrix [N_inc x nP]
% - b_theta: matrix [N_inc x nP]
% - b_phi: matrix [N_inc x nP]
% - m: vector [1 x P]
% - n: vector [1 x P]
% where nP = nMax*(nMax + 2), N_inc is the number of incidence directions
%
% Dependency:
%

theta = theta(:);
phi = phi(:);

% prefactor stuff
% ----------------
% NOTE: should eventually store those once and for all
pf = mie_prefactors(N_nmax);

% angular stuff
% ----------------
pitaurho = mie_PiTauRho(N_nmax, theta);

mphi = bsxfun(@times, pf.m, phi);
expimphi = exp(-1i*mphi);

% combining pieces
% ----------------

% temporary storage of intermediate matrices
% NOTE: the prefactor combs could be pre-calculated in spacPrefactors
% pf.in = (1i).^n;
% pf.m1 = (-1).^m;
% pf.rn = sqrt((2*n+1) ./ ( 4*pi*n.*(n+1)))
pref = pf.m1 .* pf.in .* pf.rn; % (-1)^m * i^n * r_n

%a_theta = (-1)^(m+1) i^(n+1) r_n exp(-imphi) pi_mn
stMieab.a_theta = bsxfun(@times, expimphi .* pitaurho.pinm, -1i * pref);
stMieab.a_phi   = bsxfun(@times, expimphi .* pitaurho.taunm,-1  * pref);

stMieab.b_theta =  1i * stMieab.a_phi;
stMieab.b_phi   = -1i * stMieab.a_theta;

% will need those so pass along
stMieab.m =  pf.m;
stMieab.n =  pf.n;
stMieab.nMax =  N_nmax;

end
