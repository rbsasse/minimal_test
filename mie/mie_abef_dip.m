function stMieabef = mie_abef_dip(kn, r, theta, phi, P, Ep0, stMieabang)
%% mie_abef_dip
% Calculates Mie coefficients a,b,e,f for dipole sources
%
%	mie_abef_dip(k, r, p, stMieabang) computes amn, bmn, emn, fmn
%
% Input:
% - k: scalar integer wavenumber
% - Ep0: complex scalar prefactor ik^3/(eps1*eps0) [note: no |p| (elsewhere)]
% - r: numeric vector [nR x 1] dipole positions
% - theta: numeric vector [nR x 1] dipole positions
% - phi: numeric vector [nR x 1] dipole positions
% - P: complex matrix [3Nr x Na] dipole moments
% - stMieabang: structure with geom parts of amn's bmn's for all positions
%
% Output: structure with fields:
% - amn: matrix [N_inc x nP]
% - bmn: matrix [N_inc x nP]
% - emn: matrix [N_inc x nP]
% - fmn: matrix [N_inc x nP]
%
% where N_inc is the number of columns of P
% N_dip is the number of dipoles
% nP = nMax*(nMax + 2)
%
% Dependency:
% cartesian_to_spherical, mie_Zn


%% wavelength-dependent stuff
% (NOTE we work one lambda at a time)
% maybe should calculate all at once instead
% ----------------
kr = kn*r;
N_nmax = stMieabang.nMax;
stZn = mie_Zn(kr, N_nmax);
% repeat Zn's for all m's
idn = stMieabang.n;
% Z is % [nR x nMax] -> repeat for all m's
rep_Zr0 = stZn.Zr0(:,idn);% [nR x nP]
rep_Zi0 = stZn.Zi0(:,idn);% [nR x nP]
rep_Zr1 = stZn.Zr1(:,idn);% [nR x nP]
rep_Zi1 = stZn.Zi1(:,idn);% [nR x nP]
rep_Zr2 = stZn.Zr2(:,idn);% [nR x nP]
rep_Zi2 = stZn.Zi2(:,idn);% [nR x nP]

%% dipole components
% ----------------

% loop over columns of P (corresponding to independent incident angles)
N_inc = size(P,2);
N_dip = size(P,1)/3;

for (ii=1:N_inc)

    % select one column and make it a 3xNr matrix
    P_tmp = reshape(P(:,ii),3, N_dip);

    P_norm = sqrt(sum(P_tmp.*conj(P_tmp)));
    % note: |p| will go with Ep0, here need p_hat (normalised)
    P_tmp = bsxfun(@times, P_tmp, 1./P_norm);
    % need to convert from cart to spherical
    [ P_sph ] = cartesian_to_spherical(theta, phi, P_tmp);


    P_r     = P_sph(1,:).';
    P_theta = P_sph(2,:).';
    P_phi   = P_sph(3,:).';


    %% combining pieces
    % ----------------

    % product for all dipoles
    % Eqs. 20, 21
    % a_theta is [N_dip x nP] (dipoles x m&n's)
    % p_theta is [N_dip x 1]
    % rep_Z  is [N_dip x nP]
    a_t = bsxfun(@times, stMieabang.a_theta, P_theta); % all [N_dip x nP]
    a_p = bsxfun(@times, stMieabang.a_phi, P_phi);
    a_mns =  rep_Zi0 .* (a_t + a_p);                  % irregular
    e_mns =  rep_Zr0 .* (a_t + a_p);                  % regular

    b_t = bsxfun(@times, stMieabang.b_theta, P_theta); % all [N_dip x nP]
    b_p = bsxfun(@times, stMieabang.b_phi,   P_phi);
    b_r = bsxfun(@times, stMieabang.b_r,     P_r);
    b_mns =  rep_Zi2 .* (b_t + b_p) + rep_Zi1 .* b_r; % irregular
    f_mns =  rep_Zr2 .* (b_t + b_p) + rep_Zr1 .* b_r;  % regular

    %% finally sum over dipoles (! not m,n's) to get combined excitation
    % ----------------
    Ep0_p = Ep0 * P_norm.'; % reintroduce |p|
    % ab
    stMieabef.amn(ii,:) = sum(bsxfun(@times, Ep0_p, a_mns), 1); % sum over rows (i.e. dipoles)
    stMieabef.bmn(ii,:) = sum(bsxfun(@times, Ep0_p, b_mns), 1); % sum over rows (i.e. dipoles)
    % ef
    stMieabef.emn(ii,:) = sum(bsxfun(@times, Ep0_p, e_mns), 1); % sum over rows (i.e. dipoles)
    stMieabef.fmn(ii,:) = sum(bsxfun(@times, Ep0_p, f_mns), 1); % sum over rows (i.e. dipoles)

end

stMieabef.n = stMieabang.n; % will need those so pass along
stMieabef.m = stMieabang.m;

end
