function Csca = mie_Csca_bare(stMiepq, kn)
%% mie_Csca
% Calculates the scattering cross-section

Csca = 1/kn^2*sum(abs(stMiepq.pmn).^2 + ...
    abs(stMiepq.qmn).^2, 2);


end
