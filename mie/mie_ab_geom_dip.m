function stMieab = mie_ab_geom_dip(theta, phi, N_nmax)
%% mie_ab_geom_dip
% Calculates "geometrical" part of Mie coefficients a and b for given
% dipole positions
%
%	mie_ab_geom_dip(theta, phi, nMax) computes angular amn, bmn
%
% These coefficients need to be combined with the spherical Bessel and
% spherical components of P to produce the actual a and b of the multipole expansion.
%
% Input:
% - theta: numeric vector [N_dip x 1] dipole position angle
% - phi: numeric vector [N_dip x 1] dipole position angle
% - nMax: scalar integer [1 x 1]
%
% Output: structure with fields:
% - a_theta: matrix [N_dip x nP]
% - a_phi: matrix [N_dip x nP]
% - b_theta: matrix [N_dip x nP]
% - b_phi: matrix [N_dip x nP]
% - b_r: matrix [N_dip x nP]
% - m: vector [1 x nP]
% - n: vector [1 x nP]
% where nP = nMax*(nMax + 2), N_dip is the number of dipoles
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
% pf.m1 = (-1).^m;
% pf.nnp1 = n.*(n+1);
% pf.rn = sqrt((2*n+1) ./ ( 4*pi*n.*(n+1)))
pref = pf.m1 .* pf.rn; % (-1)^m * r_n

%a_theta = (-1)^(m+1) i r_n exp(-imphi) pi_mn
stMieab.a_theta = bsxfun(@times, expimphi .* pitaurho.pinm, -1i * pref);
%a_phi = (-1)^(m+1) r_n exp(-imphi) tau_mn
stMieab.a_phi   = bsxfun(@times, expimphi .* pitaurho.taunm,-1 * pref);

stMieab.b_theta = -1 * stMieab.a_phi;
stMieab.b_phi   =      stMieab.a_theta;
%b_r = n(n+1) (-1)^m r_n exp(-imphi) rho_mn
stMieab.b_r     = bsxfun(@times, expimphi .* pitaurho.rhonm, pf.nnp1 .* pref);
% will need those so pass along
stMieab.m =  pf.m;
stMieab.n =  pf.n;
stMieab.nMax =  N_nmax;

end
