function stZn = mie_Zn(rho, N_nmax)
% computes the three Zn(rho) auxiliary functions for the radial
% dependence of VSHs for n=1 to nNmax
% used for both regular VSHs (j) and irregular VSHs (h1).
%
% Parameters:
% - N_nmax: scalar integer
%          number of n in series
% - rho:   column vector [R x 1] (no zero components allowed, even for regular VSH, for speed optimization)
%          arguments of the VSHs
%
%
% Returns: stZn structure with 6 fields
%          containing matrices [R x nNmax]
% - stZn.Zr0 is regular Z_n^0(rho)
% - stZn.Zr1 is regular Z_n^1(rho)
% - stZn.Zr2 is regular Z_n^2(rho)
% - stZn.Zi0 is irregular Z_n^0(rho)
% - stZn.Zi1 is irregular Z_n^1(rho)
% - stZn.Zi2 is irregular Z_n^2(rho)
%

if ~isempty(find(rho==0,1))
    disp 'Warning: rho=0 arguments not allowed in mie_Zn...'
end

n=1:N_nmax;
nm1=0:N_nmax;
nu=nm1+0.5;

[numat, rhomat] = meshgrid(nu,rho);

fr=besselj(numat,rhomat);
fi=besselh(numat, rhomat);


% fr and fi are matrices [R x nNmax+1] of cylindrical Bessel
% Z_{n+0.5}(rho), n=0..nNmax

sq=sqrt((pi/2)./rho); % [R x 1]
fr=bsxfun(@times,fr,sq); % [R x nNmax+1]
fi=bsxfun(@times,fi,sq); % [R x nNmax+1]
% fr and fi are now matrices of spherical Bessel
% z_n(rho), n=0..nNmax or equivalently z_{n-1}(rho), n=1..nNmax+1

stZn.Zr0=fr(:,2:(N_nmax+1));
stZn.Zi0=fi(:,2:(N_nmax+1));
stZn.Zr1=bsxfun(@times,stZn.Zr0, 1./ rho);
stZn.Zi1=bsxfun(@times,stZn.Zi0, 1./ rho);

% Computes: Z2_n=z_{n-1} - n Z1_n
stZn.Zr2 = fr(:,n) - bsxfun(@times,n,stZn.Zr1);
stZn.Zi2 = fi(:,n) - bsxfun(@times,n,stZn.Zi1);

end
