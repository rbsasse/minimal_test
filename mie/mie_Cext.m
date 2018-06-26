function [Cext] = mie_Cext(stMiepq, stMieab, stMieef, kn)
%% mie_Cext
% Calculates the extinction cross-section of the combined system using Mie
% coefficients
%
%	mie_Cext(stMiepq, stMieab) computes Cext = sum(Re(a.p* + a.e*) + Re(b.q* + b.f*))
%
% Input:
% - stMiepq: structure with p's and q's describing the sphere's total scattered field
% - stMieabef: structure with a's and b's describing incident PWE
% - stMieef: structure with e's and f's describing outgoing radiation from dipoles
% - kn: wavenumber in medium (scalar)
%
% Output: Extinction cross-section
%
% Dependency: none
%

% note: the incident field may have multiple directions, proceed one at a time
%Na = size(stMieab.amn, 1);

% Cext = zeros(1, Na);
% for (ii=1:Na)
    % select row for one incident angle
%     PWEamn = stMieab.amn(ii,:);
%     PWEbmn = stMieab.bmn(ii,:);

    % it would be nice to return contributions separately but it's not
    % directly meaningful if p and q contain the contribution from dipoles
    % scattered by the sphere... Need p,q just from PWE a,b x Gamma, Delta

    %     Cext_sphere(ii) = sum(real(PWEamn .* conj(stMiepq.pmn)) + ...
    %                            real(PWEbmn .* conj(stMiepq.qmn)));
    %     Cext_dip(ii) = sum(real(PWEamn .* conj(stMieef.emn)) + ...
    %                            real(PWEbmn .* conj(stMieef.fmn)));

    Cext = -1/kn^2 * sum(real(stMieab.amn .* conj(stMiepq.pmn + stMieef.emn)) + ...
                   real(stMieab.bmn .* conj(stMiepq.qmn + stMieef.fmn)),2);
% end


end
