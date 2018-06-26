function stRBxi2 = mie_RBxi2(N_nmax, rho)
% Calculates the Ricatti-Bessel function xi_n(rho)=h1_n(rho)*rho and its derivative for n=1:nNmax
%
% Parameters:
% - N_nmax: scalar integer
%          number of n in series
% - rho:   column vector [R x 1] (no zero components allowed)
%          arguments of the Ricatti-Bessel function
%
% Returns: structure with 2 fields
%          each a matrix [R x N_nmax]
%          with values for each rho and n=1..N_nmax
% stRBxi2.xi:  xi_n(rho)
% stRBxi2.Dxi: xi'_n(rho)}
%


n=1:N_nmax;
nm1=0:N_nmax;
nu=nm1+0.5;

[numat, rhomat] = meshgrid(nu,rho);
f=besselh(numat,rhomat);
% f is matrix [R x nNmax+1] of cylindrical Hankel
% H1_{n+0.5}(rho), n=0..nNmax

sq=sqrt((pi/2).*rho); % [R x 1]
f=bsxfun(@times,f,sq); % [R x nNmax+1]
% f is now matrix of Ricatti-Bessel
% xi_n(rho), n=0..nNmax or equivalently xi_{n-1}(rho), n=1..nNmax+1

stRBxi2.xi=f(:,n+1);

% Computes: xi_n'=xi_{n-1} - n xi_n/rho
% and check for loss of precision in sum
% Note that (1./rho) * n is a matrix [R x nMax]
stRBxi2.Dxi = f(:,n) - ( (1./rho)*n ).*f(:,n+1);

end