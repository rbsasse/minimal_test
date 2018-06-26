function [nodes, weights] = quadrature_sphere(N, method)
%QUADRATURE_SPHERE 2D quadrature over all incident directions
%
% Returns quadrature points and nodes to perform full angular averaging over incident directions
% integral(f) = weights * f(alpha, beta)'
% [nodes, weights] = quadrature_sphere(10, 'qmc')
% NOTE: 'gl' and 'fibonacci' methods may return (slightly) more than N points;
% 'cheap' always returns 3 points
%
% PARAMETERS:
% - N
% - method: type of quadrature (qmc, random, gl, fibonacci, cheap)
%
% RETURNS: quadrature points and weights
%
% DEPENDS: quadrature_lgwt
%
% FAMILY: low_level, utility, quadrature
%

switch method

case 'none' % just z axis
  nodes = [0; pi/2; 0]; % along z
  weights = [1];
  
case 'cheap' % few angles
  nodes = [0, pi/2, 0; % along z
      pi/2, 0, 0; % along x
      pi/2, pi/2, 0]'; % along y
  weights = [1 1 1]'./3;

    case 'qmc' % quasi monte-carlo
        p = haltonset(2, 'Skip', 1, 'Leap',0);
        x =  net(p, N);
        alpha = 2*pi*x(:,1);
        beta = acos(2*x(:,2) - 1); % cos(beta) in [-1,1]
        nodes = [alpha, beta, 0*beta]';
        weights = repmat(1/N, N, 1);

    case 'random' % monte-carlo with random points
        phi = pi*(2*rand(N,1) - 1); % uniform [-pi,pi]
        theta = acos(2*rand(N,1) - 1); % cos-uniform [-1,1]
        nodes = [phi, theta, 0*theta]';
        weights = repmat(1/N, N, 1);

    case 'gl' % gauss legendre
        % might have slightly more than N total points
        rndN = ceil(sqrt(N/2));
        [alphan, alphaw] = quadrature_lgwt(2*rndN, 0, 2*pi);
        [betan, betaw] = quadrature_lgwt(rndN, 0, 1);
        betan = acos(2*betan - 1); % cos(beta) in [-1,1]
        [nalpha, nbeta] = meshgrid(alphan, betan);
        [walpha, wbeta] = meshgrid(alphaw, betaw);
        nodes = [nalpha(:), nbeta(:), 0*nbeta(:)]';
        weights = 1/(2*pi)*walpha(:).*wbeta(:);

    case 'fibonacci' % fibonacci points
        No = N;
        if(mod(N, 2) == 0)
            No = N+1;
        end
        P = (No-1)/2;
        ii = (-P):P;
        colat = acos(2*ii'/No); % colatitude
        Phi = (1+sqrt(5))/2; %golden ratio
        long = 2*pi*ii' / Phi;
        nodes = [long(:), colat(:), 0*long(:)]';
        weights = repmat(1/No, No, 1);


    otherwise
        warning('Unexpected quadrature type.')
end


end
