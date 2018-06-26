function wavelength_extinction_out = compare_wavelengths_extinction(min_wavelength,step,max_wavelength,test_wavelength,n_particles)
tic;
stepz = .1;
wavelength = min_wavelength:step:max_wavelength;
epsilon= epsilon_Ag(wavelength);
material.wavelength = wavelength;
material.epsilon = epsilon;
material.medium=1.33;
ext_amp = zeros(round(pi/step)+1,2);
count = 1;
max = 0;
out = [];

for dh = 0:stepz:pi  %go through each angle orientation
%since only changing one and objects are symmetrical only need up to pi

% cluster
d = 100; %distance between nanoparticles in nanomemters 
a = 35; b=12; c=12; %set the size of the nanoparticles, ellipsoid, same sizes
dihedral = dh; alpha1 = 0*pi/180; alpha2 = pi/6;%0;
cl = cluster_dimer(d, a,b,c, dihedral, alpha1, alpha2); %orient the particles 

% %Cluster V2 adding more particles
% a=35; b=12; c=12;
% % dimer
% cl.positions = [0 0 0; 0 80 0; 80 0 0]'; % particle positions %added third row
% cl.angles    = [dh 0 0;   0 0 0; 0 0 0]'; % particle orientations
% cl.sizes     = [a b c;   a b c; a b c]'; % particle sizes

%disp('***');
xsec = spectrum_oa(cl, material);
%disp('***');
activity_extinction = xsec(:,5);
ext_amp(count,1) = dh;
ext_amp(count,2) = activity_extinction(round((test_wavelength-min_wavelength)/step));  %find value of the function at a selected wavelength
count = count + 1; 
end

for idx = 1:length(ext_amp(:,2)) %go through the list of areas and find the max
    if ext_amp(idx,2) > max
        max = ext_amp(idx,2);
        out = [ext_amp(idx,1),max];
    end
end
plot(ext_amp(:,1),ext_amp(:,2)); %plot the different amplitudes vs angles to see optimization
xlabel('Angle (radians)');
ylabel('Value of Extinction Curve (nm^2)');
title(['Dihedral Angle Orientation vs Optical Activity Extinction at ' num2str(test_wavelength) 'nm']);
toc;
wavelength_extinction_out = out;
