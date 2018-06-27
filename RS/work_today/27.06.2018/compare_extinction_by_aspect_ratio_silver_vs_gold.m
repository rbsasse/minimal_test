function extinction_out = compare_extinction_by_aspect_ratio_silver_vs_gold(min_wavelength,step,max_wavelength)
%look at extinction transmission etc as function of size, aspect ratio, and
% material
% How many particles should be in this simulation, and how should they be
% oriented etc.
tic;

wavelength = min_wavelength:step:max_wavelength;
epsilon= epsilon_Ag(wavelength);
material.wavelength = wavelength;
material.epsilon = epsilon;
material.medium=1.33;




%scale_size = 10e-1; %set a scaler
upper_lim = 5;
lower_lim = 1;
step = (upper_lim-lower_lim)/10;

d = 80; %separation
a=d/upper_lim/2; b=d/upper_lim/2; c=d/upper_lim/2;


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
cl.positions = [0 0 0; 0 d 0; d 0 0]'; % particle positions %added third row
cl.angles    = [pi/4 0 0; pi/4 0 0;pi/4 0 0]'; % particle orientations
cl.sizes     = [scale*a b c; scale*a b c; scale*a b c]'; % particle sizes

%disp('***');
xsec = spectrum_oa(cl, material);

figure()
set(groot,'defaultAxesColorOrder', visual_colorscale(3), 'defaultLineLineWidth', 2)

subplot(4,1,1)
plot(wavelength, xsec(:,1:3));
xlabel('Wavelength /nm');
ylabel('Optical cross-section /nm^2');
title(['Average spectra Silver: '  num2str(scale*a/b) ':1 aspect ratio nanoparticles']);
legend({'extinction','absorption','scattering'});

subplot(4,1,2)
plot(wavelength, 2*xsec(:,5:7));
xlabel('Wavelength /nm');
ylabel('Optical cross-section /nm^2');
title(['Optical activity Silver: ' num2str(scale*a/b) ':1 aspect ratio nanoparticles']);

epsilon= epsilon_Au(wavelength);
material.epsilon = epsilon;

xsec = spectrum_oa(cl, material);

set(groot,'defaultAxesColorOrder', visual_colorscale(3), 'defaultLineLineWidth', 2)

subplot(4,1,3)
plot(wavelength, xsec(:,1:3));
xlabel('Wavelength /nm');
ylabel('Optical cross-section /nm^2');
title(['Average spectra Gold: ' num2str(scale*a/b) ':1 aspect ratio nanoparticles']);
legend({'extinction','absorption','scattering'});

subplot(4,1,4)
plot(wavelength, 2*xsec(:,5:7));
xlabel('Wavelength /nm');
ylabel('Optical cross-section /nm^2');
title(['Optical activity Silver: ' num2str(scale*a/b) ':1 aspect ratio nanoparticles']);
figure()
visual_cluster(cl);

end


toc;
extinction_out = 0;
