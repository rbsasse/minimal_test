function [E, P, niter, error] = iterate_field(E0, G, kn, AlphaBlocks, tol, maxiter)
% OOS solution
% iterate until extinction is not converging within tolerance or reached
% maxiter

error = Inf;
niter = 1;
Etmp = E0;
E = E0;
P = polarization(E0, AlphaBlocks);
cext = c_ext(kn, P, E0); 

while((error > tol) && (niter < maxiter))
    
    Etmp = G * Etmp;
    E = E + Etmp;
    P = polarization(E, AlphaBlocks);
    cextnew = c_ext(kn, P, E0);
    error = max(abs((cextnew - cext)) / abs(cextnew + cext));
    cext = cextnew;
    niter = niter+1;
    fprintf('iteration %i, error %.2e \n', niter,error);
    
end

end

