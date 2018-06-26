function wavelength_extinction_out = compare_wavelengths_extinction_choose_particle_num(min_wavelength,step,max_wavelength,test_wavelength,n_particles)
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


%Cluster V2 adding more particles
a=35; b=12; c=12;
% dimer
cl.positions = zeros(n_particles,3)'; % particle positions %added third row
cl.angles    = zeros(n_particles,3)'; % particle orientations
cl.sizes     = zeros(n_particles,3)'; % particle sizes
cl.sizes(:,1) = a;
cl.sizes(:,2) = b;
cl.sizes(:,3) = c;


if n_particles >= 6 %6 or more particles
    
for idx = 2:2: round(n_particles/3) %place particles around on the axis separated by 100 (for more than 6 particles
    cl.positions(idx,1) = idx*100;
    cl.positions(idx+1,1) = -1*idx*100;
end
for idx = round(n_particles/3):2: round(2*n_particles/3)
    cl.positions(idx,2) = idx*100;
    cl.positions(idx+1,2) = -1*idx*100;
end
for idx = round(2*n_particles/3):2: n_particles
    cl.positions(idx,3) = idx*100;
    if idx+1 <= n_particles
    cl.positions(idx+1,3) = -1*idx*100;
    end
end

else 
% if less than 6 particles position intialization
if n_particles == 2
cl.positions(2,1) = 100;
end
if n_particles == 3
    cl.positions(2,1) = 100;
    cl.positions(3,1) = -100;
end
if n_particles == 4
    cl.positions(2,1) = 100;
    cl.positions(3,1) = -100;
    cl.positions(4,2) = 100;
end
if n_particles == 5
    cl.positions(2,1) = 100;
    cl.positions(3,1) = -100;
    cl.positions(4,2) = 100;
    cl.positions(5,2) = 100;
end
end
disp(cl);
for dh = 0:stepz:pi  %go through each angle orientation
%since only changing one and objects are symmetrical only need up to pi

% % cluster
% d = 100; %distance between nanoparticles in nanomemters 
% a = 35; b=12; c=12; %set the size of the nanoparticles, ellipsoid, same sizes
% dihedral = dh; alpha1 = 0*pi/180; alpha2 = pi/6;%0;
% cl = cluster_dimer(d, a,b,c, dihedral, alpha1, alpha2); %orient the particles 

cl.angles(1,1) = dh %changining only one angle of one particle right now 
%%%%%%% Will update code this afternoon to change this code (22.06.2018)


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
