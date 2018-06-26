function wavelength_extinction_out = compare_wavelengths_extinction_4_spheres(min_wavelength,step,max_wavelength,test_wavelength)
tic;

wavelength = min_wavelength:step:max_wavelength;
epsilon= epsilon_Ag(wavelength);
material.wavelength = wavelength;
material.epsilon = epsilon;
material.medium=1.33;
ext_amp = zeros(round(pi/step)+1,2);
count = 1;
max = 0;
out = [];
abs_area = 0; %initialize

%Cluster V2 adding more particles
a=10; b=10; c=10; %spheres
% dimer
cl.positions = zeros(4,3)'; % particle positions %added third row
cl.angles    = zeros(4,3)'; % particle orientations
cl.sizes     = zeros(4,3)'; % particle sizes
cl.sizes(:,1) = a;
cl.sizes(:,2) = b;
cl.sizes(:,3) = c;
cl.angles    = zeros(4,3)'; % particle orientations
param_min = -30;
param_max = 30
param_step = 10
disp(cl);



for  x1pos = param_min:param_step:param_max %go through positioning x,y and z particle 1
    cl.positions(1,1) = x1pos;
    for y1pos = param_min:param_step:param_max 
        cl.positions(1,2) = y1pos;
        for z1pos = param_min:param_step:param_max %three coordinates particle 1
            cl.positions(1,3) = z1pos;
            for x2pos = param_min:param_step:param_max
                cl.positions(2,1) = x2pos;
                for y2pos = param_min:param_step:param_max
                     cl.positions(2,2) = y2pos;
                    for z2pos = param_min:param_step:param_max
                         cl.positions(2,3) = z2pos;
                        for x3pos = param_min:param_step:param_max
                             cl.positions(3,1) = x3pos;
                            for y3pos = param_min:param_step:param_max
                                 cl.positions(3,2) = y3pos;
                                for z3pos = param_min:param_step:param_max
                                     cl.positions(3,3) = z3pos;
                                     %disp(cl.positions);
                                     if x1pos ~=0 && y1pos~= 0 && z1pos~= 0 && x2pos ~=0 && y2pos~= 0 && z2pos~= 0 && x3pos ~=0 && y3pos~= 0 && z3pos~= 0 && x1pos ~= x2pos && x1pos ~= x3pos && x2pos ~= x3pos && y1pos ~= y2pos && y1pos ~= y3pos && y2pos ~= y3pos && z1pos ~= z2pos && z1pos ~= z3pos && z2pos ~= z3pos
                                           %disp('***');
                                            disp('***');
                                            xsec = spectrum_oa(cl, material);
                                            disp(xsec);
                                            activity_extinction = xsec(:,5);
                                            abs_area = sum(abs(activity_extinction))
                                            if abs_area > max
                                                max = abs_area
                                                out = cl.positions
                                            end
                                     end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    visual_cluster(cl);
end

% % cluster
% d = 100; %distance between nanoparticles in nanomemters 
% a = 35; b=12; c=12; %set the size of the nanoparticles, ellipsoid, same sizes
% dihedral = dh; alpha1 = 0*pi/180; alpha2 = pi/6;%0;
% cl = cluster_dimer(d, a,b,c, dihedral, alpha1, alpha2); %orient the particles 




% plot(ext_amp(:,1),ext_amp(:,2)); %plot the different amplitudes vs angles to see optimization
% xlabel('Angle (radians)');
% ylabel('Value of Extinction Curve (nm^2)');
% title(['Dihedral Angle Orientation vs Optical Activity Extinction at ' num2str(test_wavelength) 'nm']);
toc;
wavelength_extinction_out = [out]
