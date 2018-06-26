function stMiepq = mie_pq(stSuscep, stMieab)
%% mie_pq
% Calculates Mie coefficients p and q for scattered field
%
%	mie_pq(stSuscep, stMieab) computes pmn, qmn
%
% Input:
% - stSuscep: structure with sphere susceptibilities Gamma, Delta at that wavelength
% - stMieab: structure with a's and b's describing complete incident field
%
% Output: structure with fields:
% - pmn: matrix [nR x nP]
% - qmn: matrix [nR x nP]
% where nP = nMax*(nMax + 2), nR
% The arrays use the p-index which stores the possible values
% of(n,m) in a linear array using the following convention p = n*(n+1)/2.
%
% Dependency:
%

% NOTE: Delta, Gamma calculated separately for all lambda at once
% stSuscep=mie_GDAB(nNmax,s,x);


% replicate Delta_n and Gamma_n for all m's

repDelta = stSuscep.Delta(stMieab.n); % index with repeats
repGamma = stSuscep.Gamma(stMieab.n); % index with repeats

%% return the scattered field coefficients
% ----------------

% stMiepq.pmn = repGamma .* stMieab.amn;
% stMiepq.qmn = repDelta .* stMieab.bmn;

stMiepq.pmn = bsxfun(@times,repGamma,stMieab.amn);
stMiepq.qmn = bsxfun(@times,repDelta,stMieab.bmn);

end
