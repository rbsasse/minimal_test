function wavelength_extinction_out = compare_wavelengths_extinction_4_spheres_V2(min_wavelength,step,max_wavelength,test_wavelength)
tic;


wavelength = min_wavelength:step:max_wavelength;
epsilon= epsilon_Ag(wavelength);
material.wavelength = wavelength;
material.epsilon = epsilon;
material.medium=1.33;
output = [];

 
%Cluster V2 adding more particles
a=5; b=5; c=5;
% dimer
cl.positions = [0 0 0; 0 0 0; 0 0 0; 0 0 0]'; % particle positions %added third row
cl.angles    = [pi/4 0 0;   0 pi/8 0; pi 0 pi/3; pi 0 pi/3]'; % particle orientations
cl.sizes     = [a b c;   a b c; a b c; a b c]'; % particle sizes
%xsec = spectrum_oa(cl, material);
ext_amp = zeros(round(pi/step)+1,2);
count = 1;
max = 0;
max_act_ext = [];
out = [];
abs_area = 0; %initialize

%Cluster V2 adding more particles
%a=10; b=10; c=10; %spheres


% dimer

param_min = -30;
param_max = 30;
param_step = 10;
disp(cl);



for  x1pos = param_min : param_step : param_max %go through positioning x,y and z particle 1
    for y1pos = param_min : param_step : param_max 
        for z1pos = param_min : param_step : param_max %three coordinates particle 1
            for x2pos = param_min : param_step : param_max
                for y2pos = param_min : param_step : param_max
                    for z2pos = param_min : param_step : param_max
                        for x3pos = param_min : param_step : param_max
                            for y3pos = param_min : param_step : param_max
                                for z3pos = param_min : param_step : param_max
                                     %disp(cl.position);
                                     if x1pos ~=0 && y1pos~= 0 && z1pos~= 0 && x2pos ~=0 && y2pos~= 0 && z2pos~= 0 && x3pos ~=0 && y3pos~= 0 && z3pos~= 0 && x1pos ~= x2pos && x1pos ~= x3pos && x2pos ~= x3pos && y1pos ~= y2pos && y1pos ~= y3pos && y2pos ~= y3pos && z1pos ~= z2pos && z1pos ~= z3pos && z2pos ~= z3pos
                                         %disp('***');
                                           %disp('***');
                                           cl.positions = [0 0 0;x1pos y1pos z1pos; x2pos y2pos z2pos; x3pos y3pos z3pos]'; % particle positions %added third row
                                            xsec = spectrum_oa(cl, material);
                                            %disp('***');
                                            activity_extinction = xsec(:,5);
                                            abs_area = sum(abs(activity_extinction));
                                            output(1,count) = abs_area;
                                            output(2,count) = average_particle_separation(cl.positions);
                                            count = count +1;
                                            if abs_area > max
                                                max = abs_area;
                                                out = cl.positions;
                                                max_act_ext = activity_extinction;
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
end
figure(1);
plot(wavelength,2*max_act_ext);
xlabel('Wavelength /nm')
ylabel('Optical cross-section /nm^2')
title('Optical activity')
figure(2);
cl.positions = out;
visual_cluster(cl);
figure(3);
plot(output(2,:),output(1,:),'b .','MarkerSize',20);
title('Optical Activity vs Particle Density');
ylabel('Area under Optical cross-section Extinction Curve /nm^3');
xlabel('Average distance between particles (nm)');
disp('***');
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
wavelength_extinction_out = max;%[out]
