function mn = mie_mn(N_nmax)
%% mie_mn
% Calculates mn for given N_nmax
%
%	mie_mn(N_nmax) computes n,m for n=1..N_nmax, |m|<=n
%
% Input:
% - N_nmax: scalar integer [1 x 1]
%
% Output: matrix with rows:
% - m: vector [1 x P]
% - n: vector [1 x P]
%
% Dependency:
% none

% not very fast... bad matlab trick
% cf http://stackoverflow.com/q/37955633/471093
% mn = cell2mat(arrayfun(@(n) [(-n:n) ;(-n:n)*0+n], 1:nMax, 'UniformOutput', false));


u = N_nmax^2+2*N_nmax;
mn = [ ones(1,u); zeros(1,u) ];
vv = transpose(1:N_nmax);
ww = vv.^2;
mn(1, ww) = -2*vv+1;
mn(2, ww) = 1;
mn = cumsum(mn, 2);

end
