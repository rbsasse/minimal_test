function [r_dip, theta_dip, phi_dip] = dipoles_to_spherical(positions)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

r_dip = sqrt(sum(positions.^2, 1)); % r == R0 + d here for all dipoles
theta_dip = real(acos(positions(3,:) ./ r_dip)); % make sure no complex type
phi_dip   = atan2(positions(2,:), positions(1,:));
% make sure always column vectors
r_dip = r_dip(:);
theta_dip = theta_dip(:);
phi_dip = phi_dip(:);
end

