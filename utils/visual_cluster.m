function visual_cluster( cluster, colour )
%VISUAL_CLUSTER 3D plot of a cluster to help visualisation
%
% Draws particles as ellipsoids in a 3D scene with arbitrary positions,
% sizes and orientations
%
% PARAMETERS:
% - cluster
% - colour [optional] rgb matrix
%
% RETURNS: side-effect
%
% DEPENDS: visual_ellipsoids
%
% FAMILY: user_level, graphics
%

if(nargin <2)
colour = [0 0 1];
end
   X = cluster.positions(1,:);
   N_dip = length(X);
if (size(colour, 1) < N_dip)
    colour = ones(N_dip,1)*colour(1,:);
end

% ensure no size is 0, otherwise invisible
aa = cluster.sizes(1,:);
bb = cluster.sizes(2,:);
cc = cluster.sizes(3,:);
maxi = max([aa, bb, cc]);
mini = maxi/10;
aa = max(aa, mini);
bb = max(bb, mini);
cc = max(cc, mini);

visual_ellipsoids(cluster.positions(1,:),cluster.positions(2,:),cluster.positions(3,:), ...
    aa, bb, cc, ...
    cluster.angles(1,:),cluster.angles(2,:),cluster.angles(3,:), colour)
axis vis3d
axis off
rotate3d on
% 
% light('Position',[1 1 1],'Style','infinit','Color',[1 1 1]);
% light('Position',-1*[1 1 1],'Style','infinit','Color',[1 1 1]);
% lighting gouraud;

end
