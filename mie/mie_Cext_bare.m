function [Cext] = mie_Cext_bare(stMiepq, stMieab, kn)
%% mie_Cext
% Calculates the extinction cross-section

    Cext = -1/kn^2 * sum(real(stMieab.amn .* conj(stMiepq.pmn)) + ...
                   real(stMieab.bmn .* conj(stMiepq.qmn)),2);


end
