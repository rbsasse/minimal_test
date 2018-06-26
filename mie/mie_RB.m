function stRB = mie_RB(N_nmax, rho)
% Calculates both Ricatti-Bessel functions psi_n(rho) and xi_n(rho) and their derivative for n=1:nNmax
% This function is faster than calling RBpsi2 and RBxi2 sequentially
%
% Parameters:
% - N_nmax: scalar integer
%          number of n in series
% - rho:   column vector [R x 1] (no zero components allowed)
%          arguments of the Ricatti-Bessel function
%
% Returns: structure stRB with 4 fields
%             each a matrix [R x nNmax] for each rho and n=1..nNmax
% - stRB.psi:  psi_n(rho)
% - stRB.xi:   xi_n(rho)
% - stRB.Dpsi: psi'_n(rho)
% - stRB.Dxi:  xi'_n(rho)
%


n=1:N_nmax;
nm1=0:N_nmax;
nu=nm1+0.5;

[numat, rhomat] = meshgrid(nu,rho);
fj=besselj(numat,rhomat);
f=besselh(numat,rhomat);

sq=sqrt((pi/2).*rho); % [R x 1]
fj=bsxfun(@times,fj,sq); % [R x nNmax+1]
% fj is now matrix of spherical Ricati-Bessel
% psi_n(rho), n=0..nNmax or equivalently psi_{n-1}(rho), n=1..nNmax+1
f=bsxfun(@times,f,sq); % [R x nNmax+1]
% f is now matrix of Ricatti-Bessel
% xi_n(rho), n=0..nNmax or equivalently xi_{n-1}(rho), n=1..nNmax+1

% Computes: psi_n'=psi_{n-1} - n psi_n/rho
% and check for loss of precision in sum
% Note that (1./rho) * n is a matrix [R x nMax]
stRB.Dpsi=fj(:,n) - ( (1./rho)*n ).*fj(:,n+1);

% Computes: xi_n'=xi_{n-1} - n xi_n/rho
% and check for loss of precision in sum
stRB.Dxi=f(:,n) - ( (1./rho)*n ).*f(:,n+1);

stRB.psi=fj(:,n+1);
stRB.xi=f(:,n+1);

end