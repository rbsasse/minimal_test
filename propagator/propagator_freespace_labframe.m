function A = propagator_freespace_labframe(R, kn, AlphaBlocks)
%PROPAGATOR_FREESPACE_LABFRAME Fast version of the interaction matrix for free-space coupling
%
% Free-space interaction matrix, operating in the global lab frame and in cartesian coordinates.
% This corresponds to A = I - G alpha where G is the free-space Green's
% function.
% This fast version uses aggressive Matlab optimisation, and is quite unreadable.
% See the slow version to relate to what the code actually does.
%
% PARAMETERS:
% - R
% - kn: wavevector in medium (i.e. = k0*sqrt(eps_m))
% - AlphaBlocks
%
% RETURNS: 3Nrx3Nr interaction matrix for free-space coupling
%
% DEPENDS: none
%
% FAMILY: low_level, propagator
%

N_dip = size(R,2);

krjkxN = bsxfun(@minus,kn*R(1,:).',kn*R(1,:)); % This gives [N x N]
krjkyN = bsxfun(@minus,kn*R(2,:).',kn*R(2,:)); % This gives [N x N]
krjkzN = bsxfun(@minus,kn*R(3,:).',kn*R(3,:)); % This gives [N x N]
k2rjk2N = krjkxN.^2+krjkyN.^2+krjkzN.^2;
k2rjk2N(1:N_dip+1:N_dip*N_dip)=1; % To avoid problems on the diagonal
krjkN = realsqrt(k2rjk2N);
GphaseN = kn^3*exp(1i*krjkN)./krjkN;
GfactN = (1i*krjkN-1)./k2rjk2N;

njkxN = krjkxN./krjkN;
njkyN = krjkyN./krjkN;
njkzN = krjkzN./krjkN;

ncrossnNxx = njkxN.*njkxN;
ncrossnNyx = njkyN.*njkxN;
ncrossnNzx = njkzN.*njkxN;
ncrossnNyy = njkyN.*njkyN;
ncrossnNzy = njkzN.*njkyN;
ncrossnNzz = njkzN.*njkzN;

GphaseN(1:N_dip+1:N_dip*N_dip)=0; % Enforce no self-reaction
GIFact=-GphaseN.* (1+GfactN);
GcrFact=GphaseN.*(1+3*GfactN);

% G3Nxx = GcrFact.*ncrossn3Nxx + GIFact;
% G3Nyx = GcrFact.*ncrossn3Nyx;
% G3Nzx = GcrFact.*ncrossn3Nzx;
% etc...
%

ind=0:3:(3*N_dip-1);
A = eye(3*N_dip);
A(1+ind,1+ind) = A(1+ind,1+ind) + GcrFact .* (bsxfun(@times,ncrossnNxx,AlphaBlocks(1+ind,1).') + ...
    bsxfun(@times,ncrossnNyx,AlphaBlocks(2+ind,1).') +  bsxfun(@times,ncrossnNzx,AlphaBlocks(3+ind,1).')) ...
    + bsxfun(@times,GIFact,AlphaBlocks(1+ind,1).'); % Axx

A(2+ind,1+ind) = A(2+ind,1+ind) + GcrFact .* (bsxfun(@times,ncrossnNyx,AlphaBlocks(1+ind,1).') + ...
    bsxfun(@times,ncrossnNyy,AlphaBlocks(2+ind,1).') +  bsxfun(@times,ncrossnNzy,AlphaBlocks(3+ind,1).')) ...
    + bsxfun(@times,GIFact,AlphaBlocks(2+ind,1).'); % Ayx

A(3+ind,1+ind) = A(3+ind,1+ind) + GcrFact .* (bsxfun(@times,ncrossnNzx,AlphaBlocks(1+ind,1).') + ...
    bsxfun(@times,ncrossnNzy,AlphaBlocks(2+ind,1).') +  bsxfun(@times,ncrossnNzz,AlphaBlocks(3+ind,1).')) ...
    + bsxfun(@times,GIFact,AlphaBlocks(3+ind,1).'); % Azx

A(2+ind,2+ind) = A(2+ind,2+ind) + GcrFact .* (bsxfun(@times,ncrossnNyx,AlphaBlocks(1+ind,2).') + ...
    bsxfun(@times,ncrossnNyy,AlphaBlocks(2+ind,2).') +  bsxfun(@times,ncrossnNzy,AlphaBlocks(3+ind,2).')) ...
    + bsxfun(@times,GIFact,AlphaBlocks(2+ind,2).'); % Ayy

A(3+ind,2+ind) = A(3+ind,2+ind) + GcrFact .* (bsxfun(@times,ncrossnNzx,AlphaBlocks(1+ind,2).') + ...
    bsxfun(@times,ncrossnNzy,AlphaBlocks(2+ind,2).') +  bsxfun(@times,ncrossnNzz,AlphaBlocks(3+ind,2).')) ...
    + bsxfun(@times,GIFact,AlphaBlocks(3+ind,2).'); % Azy

A(3+ind,3+ind) = A(3+ind,3+ind) + GcrFact .* (bsxfun(@times,ncrossnNzx,AlphaBlocks(1+ind,3).') + ...
    bsxfun(@times,ncrossnNzy,AlphaBlocks(2+ind,3).') +  bsxfun(@times,ncrossnNzz,AlphaBlocks(3+ind,3).')) ...
    + bsxfun(@times,GIFact,AlphaBlocks(3+ind,3).'); % Azz

A(1+ind,2+ind) = A(1+ind,2+ind) + GcrFact .* (bsxfun(@times,ncrossnNxx,AlphaBlocks(1+ind,2).') + ...
    bsxfun(@times,ncrossnNyx,AlphaBlocks(2+ind,2).') +  bsxfun(@times,ncrossnNzx,AlphaBlocks(3+ind,2).')) ...
    + bsxfun(@times,GIFact,AlphaBlocks(1+ind,2).'); % Axy

A(1+ind,3+ind) = A(1+ind,3+ind) + GcrFact .* (bsxfun(@times,ncrossnNxx,AlphaBlocks(1+ind,3).') + ...
    bsxfun(@times,ncrossnNyx,AlphaBlocks(2+ind,3).') +  bsxfun(@times,ncrossnNzx,AlphaBlocks(3+ind,3).')) ...
    + bsxfun(@times,GIFact,AlphaBlocks(1+ind,3).'); % Axz

A(2+ind,3+ind) = A(2+ind,3+ind) + GcrFact .* (bsxfun(@times,ncrossnNyx,AlphaBlocks(1+ind,3).') + ...
    bsxfun(@times,ncrossnNyy,AlphaBlocks(2+ind,3).') +  bsxfun(@times,ncrossnNzy,AlphaBlocks(3+ind,3).')) ...
    + bsxfun(@times,GIFact,AlphaBlocks(2+ind,3).'); % Ayz


end
