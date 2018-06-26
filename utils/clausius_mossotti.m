function alpha = clausius_mossotti(epsilon)
  %clausius_mossotti Polarizability corresponding to a dielectric function
  %
  % Simplest form of Clausius-Mossotti.
  %
  % PARAMETERS:
  % - epsilon: dielectric function [vector]
  %
  % RETURNS: vector
  %
  % DEPENDS: none
  %
  % FAMILY: low_level, polarizability
  %

 alpha =  (epsilon - 1) ./ (epsilon + 2);

end
