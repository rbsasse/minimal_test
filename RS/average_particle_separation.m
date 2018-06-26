function ret_ave = average_particle_separation(cl_positions)
%% This function is designed to take a set of particle positions and return the average distance between particles
%% we assume that the first particles coordinates will be [0 0 0] in threespace.
%% this function can be used for any nubmer of particles but was created for 4 particles
%% The input shouls have 3 rows per column, each column gives coordinates of a differnt particle

distances = zeros(1,(factorial(length(cl_positions(1,:)))/(2*factorial(length(cl_positions(1,:))-2)))); %Binomial coefficient (n choose 2) -- number of particles factorial gives all distance combinations
%disp(factorial(length(cl_positions(1,:))))
pos = length(cl_positions(1,:)); %starting coordinate after first for loop
for idx = 2:length(cl_positions(1,:)) %distances from all particles to origin (except particle at the origin)
    distances(idx-1) = norm(cl_positions(:,idx));
end
pos = length(cl_positions(1,:)); %starting coordinate after first for loop
for num = 2:length(cl_positions(1,:))
    for ber = num+1:length(cl_positions(1,:))
        if num ~= length(cl_positions(1,:))
            distances(pos) = norm(cl_positions(:,num) - cl_positions(:,ber));
            if num ~= length(cl_positions(1,:)) -1
                pos = pos + 1; %update the position indext since there will be particle number factorial distance
            end
        end
    end
end
%disp(distances);
ret_ave = mean(distances);