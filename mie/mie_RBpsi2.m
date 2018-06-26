function stRBpsi2 = mie_RBpsi2(N_nmax, rho)
% Calculates the Ricatti-Bessel function psi_n(rho)=j_n(rho)*rho and its derivative for n=1:nNmax
%
% Parameters:
% - N_nmax: scalar integer
%          number of n in series
% - rho:   column vector [R x 1] (with possibly zero components)
%          arguments of the Ricatti-Bessel function
%
% Returns: structure with 2 fields
%          each a matrix [R x N_nmax]
%          with values for each rho and n=1..N_nmax
% stRBpsi2.psi:  psi_n(rho)
% stRBpsi2.Dpsi: psi'_n(rho)}
%

n=1:N_nmax;
nm1=0:N_nmax;
nu=nm1+0.5;
  
[numat, rhomat] = meshgrid(nu,rho);
f = besselj(numat,rhomat);
% f is matrix [R x nNmax+1] of cylindrical Bessel
% J_{n+0.5}(rho), n=0..nNmax

% find and exclude indices where rho=0
ind0 = find(rho==0); % ind0 is possibly empty vector
% replace with 1 to avoid error messages
rho(ind0) = 1;

sq=sqrt((pi/2).*rho); % [R x 1]
f=bsxfun(@times,f,sq); % [R x nNmax+1]
% f is now matrix of spherical Ricati-Bessel
% psi_n(rho), n=0..nNmax or equivalently psi_{n-1}(rho), n=1..nNmax+1

stRBpsi2.psi=f(:,n+1);

% Computes: psi_n'=psi_{n-1} - n psi_n/rho
% and check for loss of precision in sum
% Note that (1./rho) * n is a matrix [R x nMax]
stRBpsi2.Dpsi=f(:,n) - ( (1./rho)*n ).*f(:,n+1);

% Replace proper value for cases where rho=0, i.e. psi_n'(0)=0
if (~isempty(ind0))
    stRBpsi2.Dpsi(ind0,:)=zeros(length(ind0),N_nmax);
end

end

