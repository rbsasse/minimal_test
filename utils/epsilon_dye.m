function epsDye = epsilon_dye(alpha, rho, medium)

  eps0 = 8.854e-12; % F/m
  epsM = medium.^2;
  alphaMbar = clausius_mossotti(epsM) ;
  prefactor = rho / eps0 / 3;
  alphaDbar =  prefactor * (alpha + 0 + 0) / 3 ;% orientation average
  
  sumalpha = alphaMbar + alphaDbar;
  epsDye = (1 + 2*sumalpha) ./ (1 - sumalpha);

end

