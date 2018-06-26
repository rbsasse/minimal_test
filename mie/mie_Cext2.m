function Cext = mie_Cext2(Ejones, P, kn, positions, Incidence, Egrid)
%% mie_Cext2
% Calculates the extinction cross-section of the combined system using the
% optical theorem and reciprocity
%
%	mie_Cext2(p, Incidence, kn) computes Cext = Im(E_ff . E_inc*)
%
% Input:
% - p: dipole moments previously calculated [3nR x nA] complex matrix
% - Incidence: Euler angles of incidence angles [3 x nA]
% - kn: wave number in incident medium
%
% Output: Extinction cross-section
%
% Dependency: none
%

%% First, calculate Eloc for two orthogonal polarisations (via Mie)

% incidence angles
N_inc = size(Incidence, 2);
N_dip = size(positions, 2);


% Incidence corresponds to the incident field, here for the ORT we require
% the opposite direction, i.e theta -> pi - theta, phi -> phi + pi

RevIncidence = Incidence;
RevIncidence(1,:) = pi + Incidence(1,:); % phi -> phi + pi
RevIncidence(2,:) = pi - Incidence(2,:);
RevIncidence(3,:) = Incidence(3,:);
% NOTE need to check what we're doing with polarisation, no inversion

% two arbitrarily chosen orthogonal polarisations    
Evec1 = [1 0 0].';
Evec2 = [0 1 0].';
% corresponding fields at dipole positions
Efree1 = incident_field_OT(Evec1, kn, positions, Incidence);
Efree2 = incident_field_OT(Evec2, kn, positions, Incidence);
        
% fields are bound col-wise as [3nR x nA] (same as P)
Einc = incident_field_free(Ejones, kn, positions, Incidence);
Esca = incident_field_sphere(positions, Ejones, Incidence, Egrid);
Eloc1 = incident_field_sphere(positions, Evec1, RevIncidence, Egrid);
Eloc2 = incident_field_sphere(positions, Evec2, RevIncidence, Egrid);


%% Next, use reciprocity to deduce E_ff and scattering amplitude S
% E_ff =  eikR /(kappa R) S, with S = k^2 [ (p.Eloc1) e_1 + (p.Eloc2) e_2]
% and use the OT to get Cext
% Cext = 1/k * sum{Im([ (p.Eloc1) (e_1.Einc*) + (p.Eloc2) (e_2.Einc*)])}

% matrices are 3NxNa, easier to work one angle at a time
Cext = zeros(N_inc, 1);
% work one angle (column) at a time
for (ii=1:N_inc)
    
    tmp_Ei = reshape(Einc(:,ii), 3, N_dip);  % incident fields
    tmp_Es = reshape(Esca(:,ii), 3, N_dip);  % scattered field
    tmp_P = reshape(P(:,ii), 3, N_dip);      % P
    tmp_El1 = reshape(Eloc1(:,ii), 3, N_dip);% Eloc_1
    tmp_El2 = reshape(Eloc2(:,ii), 3, N_dip);% Eloc_2
    tmp_e_1 = reshape(Efree1(:,ii), 3, N_dip); % e_1 
    tmp_e_2 = reshape(Efree2(:,ii), 3, N_dip); % e_2
    EE1 = sum(tmp_e_1 .* conj(tmp_Ei));      % e_1 . E_inc *
    EE2 = sum(tmp_e_2 .* conj(tmp_Ei));      % e_2 . E_inc *
    EP1 = sum(tmp_El1 .* tmp_P);             % Eloc_1 . P
    EP2 = sum(tmp_El2 .* tmp_P);             % Eloc_2 . P

    % summing over al dipoles
    Cext(ii) = sum(EP1.*EE1 + EP2.*EE2);

end

Cext = kn^2 * imag(Cext);

end
