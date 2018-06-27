function optmz_angle = optimize_area_extinction_distance_dihedral_map(min_wavelength,step,max_wavelength)

tic;
wavelength = min_wavelength:step:max_wavelength; %initialize conditions from baptiste
epsilon= epsilon_Ag(wavelength);
material.wavelength = wavelength;
material.epsilon = epsilon;
material.medium=1.33;


count = 1;
max = 0;
out = [];

a = 35; b=12; c=12; %set the size of the nanoparticles, ellipsoid, same sizes
dihedral = 0; alpha1 = 0*pi/180; alpha2 = pi/6;%0;
[X,Y] = meshgrid(0:.1:pi,2*b:b/2:5*b);
abs_area = zeros(19,32);
for dh = 0:.1:pi  %go through each angle orientation
%since only changing one and objects are symmetrical only need up to pi

% cluster
%d = 100; %distance between nanoparticles in nanomemters 
dihedral = dh; 
for d = 2*b:b/6:5*b
cl = cluster_dimer(d, a,b,c, dihedral, alpha1, alpha2); % map distance and angle 

%disp('***');
xsec = spectrum_oa(cl, material);
%disp('***');
activity_extinction = xsec(:,5);
%X is angle and  complete row
%Y is distance and complete columns 
abs_area(((d-(2*b))/(b/6))+1,round(dh/.1)+1) = sum(abs(activity_extinction));
%find the abslute value of area under the curve
count = count + 1; 
end
end

[X,Y] = meshgrid([0:.1:pi],[2*b:b/6:5*b]);
disp('X');
disp(size(X));
disp('Y');
disp(size(Y));
disp('abs_area')
disp(size(abs_area));
surf(X,Y,abs_area); %plot the different areas vs angles to see optimization
xlabel('Angle (radians)');
ylabel('Separation nm');
zlabel('Absolute Value of Area Under Extinction Curve (nm^3)');
title('Dihedral Angle Orientation vs Optical Activity Extinction');
toc;

optmz_angle = out;
