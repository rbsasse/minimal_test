function stMieab = mie_ab_PWE(E_sph, stMieabang)
%% mie_ab_PWE
% Calculates Mie coefficients a_mn and b_mn for PWE incident along theta, phi
%
%	mie_ab_PWE(e_sph, stMieabang) computes amn, bmn
%
% Input:
% - E_sph: complex matrix [2 x N_inc] normalised incident fields spherical components
% - stMieabang: structure with angular parts of amn's bmn's for all incident angles
%
% Output: structure with fields:
% - amn: matrix [N_inc x nP]
% - bmn: matrix [N_inc x nP]
% where nP = nMax*(nMax + 2), N_inc is the number of incidence directions
% Note: for efficiency we compute simultaneously all angles of PW incidence. 
% However, % unlike the corresponding dipole routine, those a's and b's 
% are NOT to be added; they correspond to independent simulations.
%
% Dependency:
% cartesian_to_spherical


%% field components
% ----------------

E_theta = E_sph(1,:).';
E_phi   = E_sph(2,:).';
% incident PW doesn't have a radial component


%% combining pieces
% ----------------

% products for each angles
% Eqs. 20, 21
% a_theta is [nA x nP] (inc. angles x m&n's)
% E_theta is [1 x nA]
a_t = bsxfun(@times, stMieabang.a_theta, E_theta); % [nA x nP]
a_p = bsxfun(@times, stMieabang.a_phi, E_phi);     % [nA x nP]
a_mns = a_t + a_p; % [nA x nP]

b_t = bsxfun(@times, stMieabang.b_theta, E_theta); % [nA x nP]
b_p = bsxfun(@times, stMieabang.b_phi,   E_phi);   % [nA x nP]
b_mns = b_t + b_p; % [nA x nP]

%% don't sum anything, since those are independent PWEs
% ----------------
stMieab.amn = 4*pi * a_mns; % [nA x nP]
stMieab.bmn = 4*pi * b_mns; % [nA x nP]

% will need those so pass along
stMieab.n = stMieabang.n; % [1 x nP]
stMieab.m = stMieabang.m; % [1 x nP]

end
