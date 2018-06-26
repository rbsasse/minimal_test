function [G, filter] = propagator_freespace_labframe_slow(R, kn, AlphaBlocks, dfilter)
%PROPAGATOR_FREESPACE_LABFRAME_SLOW Slow version of the interaction matrix for free-space coupling
%
% Free-space interaction matrix, operating in the global lab frame and in cartesian coordinates.
% This slow version corresponds to the direct implementation, but is slow.
% In actual code, use the fast version.
%
% PARAMETERS:
% - R
% - kn
% - AlphaBlocks
%
% RETURNS: 3Nrx3Nr interaction matrix for free-space coupling
%
% DEPENDS: none
%
% FAMILY: low_level, propagator
%

if(nargin < 4)
    dfilter = 1;
end
    

N_dip = size(R,2);
G = eye(3*N_dip);
I3 = eye(3);
    
if(nargout >1)
    return_filter = true;
    nf = 0;
    filter = zeros(9*N_dip^2 - 9*N_dip, 2); % max number of off-diag elts
else 
    return_filter = false;
end


    mind = Inf;

for jj=1:N_dip
    jjind = (jj-1)*3+1:jj*3;
    alphajj =  AlphaBlocks(jjind,1:3);
    for kk=(jj+1):N_dip
        kkind = (kk-1)*3+1:kk*3;
        alphakk =  AlphaBlocks(kkind,1:3);
        rk_to_rj = R(:,jj)-R(:,kk);
        rjk = norm(rk_to_rj);
        mind = min(mind, rjk);
        if(rjk < dfilter && return_filter)
            nf = nf+1;
            indf = (nf-1)*9+1:nf*9; % 9 elts for each pair
            [r,c] = meshgrid(jjind,kkind);
            filter(indf, :) = [r(:), c(:)];
            nf = nf+1;
            indf = (nf-1)*9+1:nf*9; % 9 elts for each pair
            filter(indf, :) = [c(:), r(:)];
        end
        rjk_hat = (rk_to_rj)/rjk;
        rjkrjk = rjk_hat*rjk_hat.';

        % note -1 prefactor cancels out as we need (I - G)
%         Gjk = exp(1i*kn*rjk)/rjk/eps_m * (kn^2*(rjkrjk - I3) + ...
%                (1i*kn*rjk-1) / rjk^2 * (3*rjkrjk - I3));

        Gjk = exp(1i*kn*rjk)/rjk * (kn^2*(rjkrjk - I3) + ...
               (1i*kn*rjk-1) / rjk^2 * (3*rjkrjk - I3));
           
        G(jjind,kkind) = Gjk   * alphakk;
        G(kkind,jjind) = Gjk.' * alphajj;
    end
end

if return_filter
    %filter(nfilter+1, :) = [mind,mind];
    %nfilter=nfilter+1;
    [r, c] = find(filter);
    filter = filter(1:max(r), :);
   %filter = filter(1:9*nf,:) ;
end

end
