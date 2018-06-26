function Csca = mie_Csca(stMiepq, stMieef, kn)
%% mie_Csca
% Calculates the scattering cross-section of the combined system
% using Mie coefficients.
%
%	mie_Csca(stMiepq, stMieef) computes Csca = sum(|e + p|^2 + |f + q|^2)
%
% Input:
% - stMiepq: structure with p's and q's describing the sphere's total scattered field
% - stMieef: structure with e's and f's describing outgoing radiation from dipoles
% - kn: wavenumber in medium (scalar)
%
% Output: Scattering cross-section
%
% Dependency: none
%


% note: rows of p,q,e,f contain independent incident directions, proceed one at a time
% Na = size(stMiepq.pmn, 1);

%for (ii=1:Na)

%     Csca(ii) = sum(abs(stMiepq.pmn(ii,:) + stMieef.emn(ii,:)).^2 + ...
%                    abs(stMiepq.qmn(ii,:) + stMieef.fmn(ii,:)).^2);

Csca = 1/kn^2*sum(abs(stMiepq.pmn + stMieef.emn).^2 + ...
    abs(stMiepq.qmn + stMieef.fmn).^2, 2);

%end

%Csca = 1/kn^2*Csca;

end
