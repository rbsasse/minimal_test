function stEAllPhi = mini_E_PWE(lambda,epsilon,cn1,dn1,r0,theta,sBessel,stPinTaun)
% Low-level function to calculate VSH expansions for PWE for a fixed r0>0, but many theta and lambda
% The VSHs used can either be of type 1 (sBessel='j') or type 3
% (sBessel='h1')
% Use r0=Inf (with 'h1') to obtain far-field properties.
% For PWE, we have: |m|=1, c_{n,1}=c_{n,-1}, d_{n,1}=-d_{n,-1}).
% The fields Ecr, Ect, Esf given in the results are discussed in the
% supplementary information.
%
% Parameters:
% - lambda: column vector [L x 1]
%           wavelengths in nm
% - epsilon: scalar or column vector [L x 1]
%           epsilon of dielectric where field is evaluated
% - cn1, dn1: 2 matrices [L x nNmax]
%             containing the Mie coefficients c_{n,1} and d_{n,1} for the
%             field expansion of M^(i)_{n,1} and N^(i)_{n,1}, respectively
% - r0:      scalar [1 x 1] possibly zero (if sBessel='j')
%           spherical coordinate r (in nm) of points
%           r0 can be Inf to obtain far-field radiation profile
% - theta:  possibly row vector [1 x T]
%           with spherical coordinate theta of points
% - sBessel: string defining the Bessel function to be used
%           sBessel='j' for or sBessel='h1'
% - stPinTaun: structure (optional)
%              with functions of theta pi_n and tau_n
%              if omitted, then the functions are computed from scrath
%              it is faster to pass this structure as argument if these
%              functions have already been calculated
%
% Returns:  stEAllPhi, structure with 3 fields
%           containing matrices [L x T]
%           representing the three components E_{cr}, E_{ct}, E_{sf} such as: of
%           E = E_{cr} cos(phi) e_r + E_{ct} cos(phi) e_theta + E_{sf} sin(phi) e_phi
% - stEAllPhi.Ecr is E_{cr}
% - stEAllPhi.Ect is E_{ct}
% - stEAllPhi.Esf is E_{sf}
%
% This file is part of the SPlaC v1.0 package (copyright 2008)
% Check the README file for further information

nNmax=size(cn1,2);
nNbLambda=length(lambda);
theta=theta(:).'; % row vector
nNbTheta=length(theta);

n=1:nNmax; %[1 x nNmax]
% need n-dep coeff for series (mu_n *n *(n+1) and mu_n)
cffnr=sqrt((2*n+1)/(4*pi)); %[1 x nNmax] for Er
mun=cffnr./(n.*(n+1)); % mu_n [1 x nNmax] for Et and Ef

if r0==0
    % special case where r0=0
    % the following is equivalent to \mathbf{E}=-4/3 * mu_1 *d11 \mathbf{e}_x
    coef1=dn1(:,1)/sqrt(3*pi); % [L x 1]
    % Results are all [L x T] matrices
    stEAllPhi.Ecr=-coef1 * sin(theta);
    stEAllPhi.Ect=-coef1 * cos(theta);
    stEAllPhi.Esf=coef1 * ones(1,nNbTheta);
    disp 'r0=0 in PweEgenThetaAllPhi';
    return
end


% get Zn(rho) for radial dependence and derived functions
if (r0==Inf)
    % for far-field radiation profile
    dn1Z1=zeros(nNbLambda,1); % [L x 1]
    icn1Z0=cn1; % [L x nNmax]
    dn1Z2=dn1; % [L x nNmax]
    mun=mun.*((-i).^(n+1));
else
    rho=2*pi* sqrt(epsilon) ./lambda * r0; % column [L x 1]
    stZnAll=mini_Z(nNmax,rho,sBessel); % fields are [L x nNmax]
    dn1Z1=dn1.*stZnAll.Z1; % [L x nNmax]
    icn1Z0=i*cn1.*stZnAll.Z0; % [L x nNmax]
    dn1Z2=dn1.*stZnAll.Z2; % [L x nNmax]
    clear stZnAll; % for memory
end
% get theta dependence if not provided
if nargin < 8
    stPinTaun=mini_pitau(nNmax,transpose(theta)); % fields are [T x nNmax]
end


% loop over lambda to avoid using 3-dimension matrices:
% At a fixed lambda, the sum over n for all theta can be carried out
% using a matrix product of the theta-and-n-dependent [T x nNmax] matrix
% by a n-dependent column vector [nNmax x 1]
% the result, a [T x 1] matrix is then transposed to a [1 x T] line
Ersum=zeros(nNbLambda,nNbTheta);
Etsum=zeros(nNbLambda,nNbTheta);
Efsum=zeros(nNbLambda,nNbTheta);
for ll=1:nNbLambda
    % for Er, vecN=d_{n,1} * Z_n^1(rho) * mu_n * n *(n+1)
    vecNdep=transpose(dn1Z1(ll,:).*cffnr); % [nNmax x 1]
    % Ersum=sum_n(pi_n(t) * vecN_n)
    % Do the sum as a matrix product (and check for loss of precision)
    Ersum(ll,:)=transpose(stPinTaun.pin * vecNdep); % [1 x T]
    % for Et and Ef
    vecNdep=transpose(icn1Z0(ll,:).*mun); % [nNmax x 1]
    vecNdep2=transpose(dn1Z2(ll,:).*mun); % [nNmax x 1]
    % Do the sums and check for loss of precision
    tmp1=stPinTaun.pin * vecNdep;
    tmp2=stPinTaun.taun * vecNdep2;
    Etsum(ll,:)=transpose(tmp1 + tmp2); % [1 x T]

    tmp1=stPinTaun.taun * vecNdep;
    tmp2=stPinTaun.pin * vecNdep2;
    Efsum(ll,:)=transpose(tmp1 + tmp2); % [1 x T]

end

% Results are all [L x T] matrices
stEAllPhi.Ecr=-2*repmat(sin(theta),nNbLambda,1).*Ersum;
stEAllPhi.Ect=-2*Etsum; % corresponds to S_2 if r0==Inf
stEAllPhi.Esf=2*Efsum; % corresponds to (-S_1) if r0==Inf

end
