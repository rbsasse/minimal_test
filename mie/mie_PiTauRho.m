function stPinmTaunm = mie_PiTauRho(N_nmax, theta)
%% mie_PiTauRho
% Calculates angular functions as defined in Mishchenko 2002
%
%	mie_PiTauRho(nMax, theta) computes angular functions pi_nm(theta)
% tau_nm(theta) and rho_nm(theta) [ Mishchenko's d^0n_0m(theta)]
% for n=1..nMax, |m|<=n, using recurrence relations.
% Those functions are defined on p. 373 of Mishchenko 2002.
%
% Input:
% - nMax: scalar integer [1 x 1]
% - theta: column vector [N_theta x 1]
%          with theta (in radians)
%          all theta's must be between 0 and pi
%
% Output: structure with fields:
% - pinm: matrix [N_theta x P] with pi_nm(theta)
% - taunm: matrix [N_theta x P] with tau_nm(theta)
% - rhomn: matrix [T x nMax] with rho_mn(theta) == d^0n_0m(theta)
% where P = nMax*(nMax + 2)
% The arrays pinm and taunm use the p-index which stores the possible values
% of(n,m) in a linear array using the following convention p = n*(n+1)/2.
%
% Dependency:
% none

if size(theta,2)>1 && size(theta,1)>1
    stop('Warning: theta must be a column vector in vshPinTaunm...');
end

theta = theta(:); % ensure column if row

if ~isempty(find(theta<0 , 1))
    disp 'Warning: theta must be >=0 in vshPinTaunm...';
end

N_theta=length(theta);
P = N_nmax*(N_nmax + 2); % maximum number of columns for pi_n and tau_n matrices
stPinmTaunm.pinm = zeros(N_theta,P);
stPinmTaunm.taunm = zeros(N_theta,P); % initialize both matrices to zero
stPinmTaunm.rhonm = zeros(N_theta,P);

% Initialize Am for m=0
Amsinmm1 = ones(N_theta,1); % [T x 1]

muc=cos(theta); % [T x 1]
mus=sin(theta); % [T x 1]

% m=0 case is treated separately
% loop on m
for m=1:N_nmax

    % Am * sin(theta)^(m-1) is computed by recurrence on m
    if m>1
        Amsinmm1 = Amsinmm1*sqrt(((2*m-1)/(2*m))).*mus;
    else
        Amsinmm1 = Amsinmm1*sqrt(((2*m-1)/(2*m)));
    end

    % Initialize recurrence pi_{m-1,m}=0, pi_{m,m}=m*Am*sin(theta)^(m-1)
    ncols = N_nmax - m + 2;
    piaux = zeros(N_theta,ncols);
    piaux(:,2) = m*Amsinmm1;
    % piaux contains pi_{m+j-2,m} j=1..(nMax-m+2),
    % i.e. pi_{m-1,m}..pi_{nMax,m}

    % Get pi_{m+1,m} to pi_nMax by recurrence
    for jj=3:ncols
        n = m + jj -2;
        % piaux(:,jj) is pi_{m+jj-2,m}
        piaux(:,jj) = (1/sqrt((n-m)*(n+m)))*((2*n-1)*muc.*piaux(:,jj-1) - sqrt((n-1-m)*(n-1+m)).*piaux(:,jj-2));
    end

    nvec = (m:N_nmax); % [1 x N2] with N2=ncols-1

    pvec = nvec.*(nvec + 1) + m;   % computes p for positive m
    pvecn = pvec - 2*m;  % computes p for negative m

    % fill in pi_nm matrix for positive or negative m
    stPinmTaunm.pinm(:,pvec) = piaux(:,2:(ncols));
    stPinmTaunm.pinm(:,pvecn)= (-1)^(m+1).*piaux(:,2:(ncols));

    % return tau_nm matrix for positive or negative m
     stPinmTaunm.taunm(:,pvec) = bsxfun(@times,piaux(:,1:(ncols-1)),-sqrt((nvec-m).*(nvec+m))/m) + (muc * (nvec/m)).*piaux(:,2:(ncols));
     stPinmTaunm.taunm(:,pvecn) = (-1)^m*stPinmTaunm.taunm(:,pvec);

    % return rho_nm matrix for positive or negative m
    stPinmTaunm.rhonm(:,pvec) = bsxfun(@times, stPinmTaunm.pinm(:,pvec), mus/m); %/m
    stPinmTaunm.rhonm(:,pvecn) = (-1)^m*stPinmTaunm.rhonm(:,pvec);

end;

% Now do m=0 case
% Initialize recurrence p_0=1, p_1=muc, t_0=0, t_1=-mus
% pnm1 contains p_{n-1} n=1..nMax+1, same for tnm1
pnm1=ones(N_theta,N_nmax+1); % [T x N+1]
pnm1(:,2)=muc; % [T x 1]
tnm1=zeros(N_theta,N_nmax+1); % [T x N+1]
tnm1(:,2)=-mus; % [T x 1]

% Get p_2 to p_nMax and t_2 to t_nMax by recurrence
% p_n is pnm1(:,n+1), t_n is tnm1(:,n+1)
for n=2:(N_nmax)
    pnm1(:,n+1)=(2*n-1)/(n)*muc.*pnm1(:,n)-(n-1)/n*pnm1(:,n-1);
    tnm1(:,n+1)=muc.*tnm1(:,n)-n*mus.*pnm1(:,n);
end;

% return p_n matrix (except n=0)
stPinmTaunm.pn0=pnm1(:,2:(N_nmax+1));
% fill in taunm for m=0 (pinm=0 for m=0 is already set)
nvec=1:N_nmax;
pvec=nvec.*(nvec+1);
stPinmTaunm.taunm(:,pvec)=tnm1(:,2:(N_nmax+1));

% Pinm = 0 for m=0, already initialized

stPinmTaunm.rhonm(:,pvec) = stPinmTaunm.pn0; % m=0
