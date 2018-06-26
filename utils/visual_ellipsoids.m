function visual_ellipsoids(X,Y,Z, a,b,c, alpha,beta,gamma, colour, nfacet)
%VISUAL_ELLIPSOIDS Makes a 3d scatter plot with ellipsoids
%
%	This function is based on SCATTER3SPH, but instead of spheres we draw ellipsoids
% which may have arbitrary orientations defined by 3 Euler angles.
%
% PARAMETERS:
% - X, Y, Z: position vectors
% - a,b,c: sizes (semi-axes of ellipsoids)
% - alpha, beta, gamma: Euler angles in ZYZ convention
% - colour [optional, Nx3 rgb matrix]
% - nfacet [optional, default 15]
% DEPENDS: none
%
% RETURNS: side-effect
%
% DEPENDS: rotation_euler_active
%
% FAMILY: high_level, graphics
%

% zero sizes don't render well, so make them 1/10 of maximum
m = max([max(a), max(b), max(c)]);
a(a == 0) = 0.1*a(a == 0);
b(b == 0) = 0.1*b(b == 0);
c(c == 0) = 0.1*c(c == 0);

% defaults
transp= 1;

if(nargin<10)
colour = ones(length(X),1)*[0 0 1];
end

if(nargin<11)
nfacet= 20;
end

%-- Sphere facets
[x, y, z] = sphere(nfacet);
original = [x(:), y(:), z(:)]';
%
[n1,n2] = size(x);

% Plot ellipsoids, in three steps:
% - scale
% - rotate
% - translate
hold on
for j= 1:length(X)

    % scale
    Scale = diag([a(j) b(j) c(j)]);
    particle = Scale * original;

    % rotate
    Rot = rotation_euler_active(alpha(j),beta(j),gamma(j));
    rs = Rot * particle; % active rotation acting on column vectors
    ssx= reshape(rs(1,:), [n2 n1])';
    ssy= reshape(rs(2,:), [n2 n1])';
    ssz= reshape(rs(3,:), [n2 n1])';

    % translate and draw
    surf(ssx+X(j), ssy+Y(j), ssz+Z(j),...
        'LineStyle','none',...
        'FaceColor',colour(j,:),...
        'Facealpha',transp);
end

daspect([1,1,1]);
