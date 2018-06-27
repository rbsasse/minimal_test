function extinction_out = compare_extinction_by_size(min_wavelength,step,max_wavelength)
%look at extinction transmission etc as function of size, aspect ratio, and
% material
tic;

wavelength = min_wavelength:step:max_wavelength;
epsilon= epsilon_Ag(wavelength);
material.wavelength = wavelength;
material.epsilon = epsilon;
material.medium=1.33;

a=35/35; b=12/35; c=12/35;

%scale_size = 10e-1; %set a scaler
step = 10;
upper_lim = 100;
lower_lim = 1;

for scale = lower_lim:step:upper_lim  %go through each possible size
epsilon= epsilon_Ag(wavelength);
material.epsilon = epsilon;


% cluster
% d = 100; %distance between nanoparticles in nanomemters 
% a = 35; b=12; c=12; %set the size of the nanoparticles, ellipsoid, same sizes
% dihedral = dh; alpha1 = 0*pi/180; alpha2 = pi/6;%0;
% cl = cluster_dimer(d, a,b,c, dihedral, alpha1, alpha2); %orient the particles 

%Cluster V2 adding more particles

% dimer
cl.positions = scale.*[0 0 0; 0 80 0; 80 0 0]'; % particle positions %added third row
cl.angles    = [0 0 0;   0 0 0; 0 0 0]'; % particle orientations
cl.sizes     = scale.*[a b c;   a b c; a b c]'; % particle sizes

%disp('***');
xsec = spectrum_oa(cl, material);

%figure()
set(groot,'defaultAxesColorOrder', visual_colorscale(3), 'defaultLineLineWidth', 2)

subplot(4,(upper_lim)/step,1+(scale-1)/step)
plot(wavelength, xsec(:,1:3));
xlabel('Wavelength /nm');
ylabel('Optical cross-section /nm^2');
title(['Average spectra Silver: ' num2str(scale) ' nm sized nanoparticle']);
legend({'extinction','absorption','scattering'});

subplot(4,(upper_lim)/step,1+(scale-1)/step+(upper_lim)/step)
plot(wavelength, 2*xsec(:,5:7));
xlabel('Wavelength /nm');
ylabel('Optical cross-section /nm^2');
title(['Optical activity Silver: ' num2str(scale) ' nm sized nanoparticles']);

epsilon= epsilon_Au(wavelength);
material.epsilon = epsilon;

xsec = spectrum_oa(cl, material);

set(groot,'defaultAxesColorOrder', visual_colorscale(3), 'defaultLineLineWidth', 2)

subplot(4,(upper_lim)/step,1+(scale-1)/step+(upper_lim)/step*2)
plot(wavelength, xsec(:,1:3));
xlabel('Wavelength /nm');
ylabel('Optical cross-section /nm^2');
title(['Average spectra Gold: ' num2str(scale) ' nm sized nanoparticles']);
legend({'extinction','absorption','scattering'});

subplot(4,(upper_lim)/step,1+(scale-1)/step +(upper_lim)/step*3)
plot(wavelength, 2*xsec(:,5:7));
xlabel('Wavelength /nm');
ylabel('Optical cross-section /nm^2');
title(['Optical activity Gold: ' num2str(scale) ' nm sized nanoparticles']);
end

toc;
extinction_out = 0;
