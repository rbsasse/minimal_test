function S = propagator_sphere_labframe(cl, wavelength, epsilon, ...
    medium, AlphaBlocks, N_esa, N_mie, N_theta, ...
    stPiTauMie, stPiTauESA)
%PROPAGATOR_SPHERE_LABFRAME Interaction matrix between dipoles mediated by a sphere
%
% The components are expressed in cartesian coordinates, in the global lab frame
%
% PARAMETERS:
% - cl cluster
% - wavelength
% - epsilon
% - medium
% - AlphaBlocks
% - nNmax
% - nMie
% - nTheta
% - stPiTauMie
% - stPiTauESA
%
% RETURNS: 3Nrx3Nr interaction matrix in lab frame
%
% DEPENDS: mie_dipole_z, rotation_euler_passive, propagator_sphere_z
%
% FAMILY: low_level, propagator, sphere
%

N_dip = size(cl.positions, 2);
S = zeros(3*N_dip, 3*N_dip); % interaction matrix
d = cl.d;
R0 = cl.R0;
R = R0+d;

if nargin<10
    theta=transpose(linspace(0,pi,N_theta)); % row [1 x T]
    stPiTauMie = mie_vshPinmTaunm01(N_mie,theta);
    stPiTauESA = mie_vshPinmTaunm01(N_esa,theta);
end


%% get the Mie solution for a fixed grid of thetas
Egrid = mie_dipole_z(d, wavelength, epsilon, medium, R0, ...
    N_esa, N_mie, N_theta, stPiTauMie, stPiTauESA);

% hpositions is a Nx3 matrix of cartesian dipole positions
% hfields is a Nx3 matrix of cartesian electric field components
% both are expressed in the rotated frame

% loop over dipoles (fill 3x3-block-columns of S)
for (jj = 1:N_dip) % source dipole
    %% rotate all positions so that dipole jj is aligned along z
    % note: this bit is purely geometrical, so could be done just once
    % for all wavelengths
    
    Rjj = sqrt(sum(cl.positions(:,jj).^2, 1)); % == R0 + d for spherical shell
    thetajj = real(acos(cl.positions(3,jj) ./ Rjj)); % acos(z/R)
    y = cl.positions(2,jj);
    x = cl.positions(1,jj);
    if(abs(x) < 1e-15 && abs(y) < 1e-15)
        phijj = 0; % no use rotating, we're very close to the North pole
    else
        phijj = atan2(y, x); % atan(y/x)
    end
    thetajj;
    phijj;
    % current 'source' dipole is oriented along phi, theta
    % to bring it along z, need passive rotation of coords by phi, theta
    %
    Rot = rotation_euler_passive(phijj, thetajj, 0);
    tRot = transpose(Rot);
    
    hpositions = Rot * cl.positions; % Rotate column vectors
    
    % rotated angular coordinates
    
    R_dip = sqrt(sum(hpositions.^2, 1)); % == R0 + d for spherical shell
    theta = real(acos(hpositions(3,:) ./ R_dip)); % make sure no complex type
    phi   = atan2(hpositions(2,:), hpositions(1,:));

    % no unnecessary rotation for the dipole itself
    phi(jj) = 0; % this value is otherwise undefined (point on z axis, atan(0/0))

    dipolejj = ((jj-1)*3+1):(jj*3);
    alphajj = AlphaBlocks(dipolejj,1:3);
    %% Green's propagator of jj to all iis, mediated by the sphere
    % hSc is S in the rotated frame where 3 sources are along z
    hSc = propagator_sphere_z(theta, phi, Egrid); % in rotated cartesian

    % rotate back and multiply by alphai
    dipolejj = ((jj-1)*3+1):(jj*3);
    for (ii = 1:N_dip) % loop over the target dipoles (including itself)
        dipoleii = ((ii-1)*3+1):(ii*3); % target, ie 3-block-row
        hSc(dipoleii, 1:3);
        Sij = tRot * hSc(dipoleii, 1:3) * Rot; % rotate back
        % SijAlphaj = Sij * alphajj
        Sij;
        S(dipoleii, dipolejj) = Sij * alphajj; % Sij and alpha are in lab frame
        %S(dipoleii, dipolejj) = Sij; % testing symmetry
    end

end


end
