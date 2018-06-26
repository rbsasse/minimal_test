function stPinmTaunm = mie_vshPinmTaunm01(N_nmax, theta)
  %mie_VSHPINMTAUNM01 Angular functions pi and tau for m=0 and m=1
  %
  % Special-case utility functionto cater for illumination along z axis
  %
  % PARAMETERS:
  % - N_nmax maximum multipole order considered
  % - theta polar angles where the functions are to be evaluated
  %
  % RETURNS: a structure with the pi and tau values at grid of positions
  %
  % DEPENDS: none
  %
  % FAMILY: low_level, mie
  %

% NOTE: only m=0 and m=1 are considered

if size(theta,2)>1
    disp 'Warning: theta must be a column vector in vshPinTaun01...';
end
if ~isempty(find(theta<0 , 1))
    disp 'Warning: theta must be >=0 in vshPinTaun01...';
end

nrows=length(theta);

muc=cos(theta); % [T x 1]
mus=sin(theta); % [T x 1]

% m=1 case
% Initialize recurrence pi_{1,0}=0, pi_{1,1}=1*A1
ncols = N_nmax +1;
piaux = zeros(nrows,ncols);
piaux(:,2) = sqrt(1/2);
% piaux contains pi_{j+1,m} j=1..(nMax+1),
% i.e. pi_{0,1}..pi_{nMax,1}

% Get pi_{2,1} to pi_nMax by recurrence
for jj=3:ncols
    n = jj - 1;
    % piaux(:,jj) is pi_{m+jj-2,m}
    piaux(:,jj) = (1/sqrt((n-1)*(n+1)))*((2*n-1)*muc.*piaux(:,jj-1) - sqrt((n-2)*n).*piaux(:,jj-2));
end;

nvec = (1:N_nmax); % [1 x N2] with N2=ncols-1

% fill in pi_n1 matrix for positive m
stPinmTaunm.pin1 = piaux(:,2:(ncols));

% return tau_n1 matrix for positive m
stPinmTaunm.taun1 = - bsxfun(@times,piaux(:,1:(ncols-1)),sqrt(nvec.^2-1)) + (muc * nvec).*piaux(:,2:(ncols));

% Now do m=0 case
% Initialize recurrence p_0=1, p_1=muc, t_0=0, t_1=-mus
% pnm1 contains p_{n-1} n=1..nMax+1, same for tnm1
pnm1=ones(nrows,N_nmax+1); % [T x N+1]
pnm1(:,2)=muc; % [T x 1]
tnm1=zeros(nrows,N_nmax+1); % [T x N+1]
tnm1(:,2)=-mus; % [T x 1]

% Get p_2 to p_nMax and t_2 to t_nMax by recurrence
% p_n is pnm1(:,n+1), t_n is tnm1(:,n+1)
for n=2:(N_nmax)
    pnm1(:,n+1)=(2*n-1)/(n)*muc.*pnm1(:,n)-(n-1)/n*pnm1(:,n-1);
    tnm1(:,n+1)=muc.*tnm1(:,n)-n*mus.*pnm1(:,n);
end;

% return p_n and tau_n0 matrices (except n=0)
stPinmTaunm.Pn0=pnm1(:,2:(N_nmax+1));
stPinmTaunm.taun0=tnm1(:,2:(N_nmax+1));

end
