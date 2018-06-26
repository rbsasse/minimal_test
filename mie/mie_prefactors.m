function stPrefactors = mie_prefactors(N_nmax)
%% mie_prefactors
% Calculates prefactors for Mie coefficients
%
%	mie_prefactors(N_nmax) computes prefactors for n=1..N_nmax, |m|<=n
%
% Input:
% - nMax: scalar integer [1 x 1]
%
% Output: structure with fields:
% - m: vector [1 x P]
% - n: vector [1 x P]
% - p: vector [1 x P]
% - m1: vector [1 x P]
% - in: vector [1 x P]
% - nnp1: vector [1 x P]
% - rn: vector [1 x P]
% where P = N_nmax*(N_nmax + 2)
% The arrays use the p-index which stores the possible values
% of(n,m) in a linear array using the following convention p = n*(n+1)/2.
%
% Dependency:
% mie_mn


mn = mie_mn(N_nmax);

m = mn(1,:);
n = mn(2,:);

stPrefactors.m = m;
stPrefactors.n = n;
stPrefactors.p = 1:size(mn,1);
stPrefactors.in = (1i).^n;
stPrefactors.m1 = (-1).^m;
stPrefactors.nnp1 = n.*(n+1);
stPrefactors.rn = sqrt((2*n+1) ./ ( 4*pi*n.*(n+1)));
% note the m=0 values are at m(nv.*(nv+1)) where nv=1:nMax

end
