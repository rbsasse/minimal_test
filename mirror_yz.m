%% This function will take in a set of coordinates of n objects and return the cordinates for the system reflected on the YZ plane
%% The input should be a matrix 3 by n 

function ret_coords = mirror_yz(positions)
mirrior_pos = positions;
mirrior_pos(1,:) = -1.*(positions(1,:));
ret_coords = mirrior_pos;
end
