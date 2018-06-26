function P = polarization(E, AlphaBlocks)
%POLARIZATION compute the polarisation P from the field E
%
% Given the cda solution E, compute P = alpha E for all dipoles and all angles of incidence
% NOTE: it could be clearer to store 3x3 tensor blocks as cube or list
% but still requires funny indexing for E and P in any case
%
% PARAMETERS:
% - E 3NrxNa matrix of incident fields
% - AlphaBlocks 3Nrx3 matrix of polarizability tensors
%
% DEPENDS: none
%
% FAMILY: low_level
%

N_dip = size(E,1)/3;
P = E; % allocate size

for jj=1:N_dip
    indjj = (jj-1)*3+1:jj*3;
    alphajj = AlphaBlocks(indjj,1:3);
    P(indjj, :) = alphajj * E(indjj, :);
end

end
