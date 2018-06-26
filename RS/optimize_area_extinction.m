function optmz_angle = optimize_area_extinction(min_wavelength,step,max_wavelength)

tic;
wavelength = min_wavelength:step:max_wavelength; %initialize conditions from baptiste
epsilon= epsilon_Ag(wavelength);
material.wavelength = wavelength;
material.epsilon = epsilon;
material.medium=1.33;
abs_area = zeros(round(pi/step)+1,2);
count = 1;
max = 0;
out = [];

for dh = 0:.01:pi  %go through each angle orientation
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
abs_area(count,1) = dh;
abs_area(count,2) = sum(abs(activity_extinction));  %find the abslute value of area under the curve
count = count + 1; 
end

for idx = 1:length(abs_area(:,2)) %go through the list of areas and find the max
    if abs_area(idx,2) > max
          max = abs_area(idx,2) ;
        out = [abs_area(idx,1),max];
    end
end
plot(abs_area(:,1),abs_area(:,2)); %plot the different areas vs angles to see optimization
xlabel('Angle (radians)');
ylabel('Absolute Value of Area Under Extinction Curve (nm^3)');
title('Dihedral Angle Orientation vs Optical Activity Extinction');
toc;

optmz_angle = out;
