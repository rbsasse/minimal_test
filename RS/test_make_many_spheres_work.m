clear
close
wavelength = 400:2:800;
epsilon= epsilon_Ag(wavelength);
material.wavelength = wavelength;
material.epsilon = epsilon;
material.medium=1.33;


% % cluster
% d = 50; %distance between nanoparticles in nanomemters 
% a = 35; b=12; c=12; %set the size of the nanoparticles, ellipsoid, same sizes
% dihedral = pi/8; alpha1 = 0*pi/180; alpha2 = 0;%0;
% cl = cluster_dimer(d, a,b,c, dihedral, alpha1, alpha2); %orient the particles 
%Cluster V2 adding more particles
a=10; b=10; c=10;
% dimer
cl.positions = [0 0 0; 0 80 0; 80 0 0; 0 0 80]'; % particle positions %added third row
cl.angles    = [pi/4 0 0;   0 pi/8 0; pi 0 pi/3; pi 0 pi/3]'; % particle orientations
cl.sizes     = [a b c;   a b c; a b c; a b c]'; % particle sizes
xsec = spectrum_oa(cl, material);



% particle positions %added third row
cl.angles    = [pi/4 0 0;   0 pi/8 0; pi 0 pi/3; pi 0 pi/3]'; % particle orientations
cl.sizes     = [a b c;   a b c; a b c; a b c]'; % particle sizes
xsec = spectrum_oa(cl, material);
