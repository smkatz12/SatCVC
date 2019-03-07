clear
close all

data = xlsread('Positions.csv');
n = 6; %number of robots
[length, width] = size(data);
timesteps = length/n; %number of timesteps
for i=1:1:timesteps
    %trajectory_angles is an n by 2*timesteps matrix containing time history of
    %all robots
    current_angles = data(i*6-5:i*6,:);
    trajectory_angles(1:n,i*2-1:i*2) = current_angles;
    for j=1:1:n
    %get trajectory_XYZ, a n by 3*timestep matrix containing time history
    %of all robots in XYZ coordinates
        phi = current_angles(j,1);
        lambda = current_angles(j,2);
        trajectory_XYZ(j,i*3-2) = cosd(phi) * cosd(lambda); %x position
        trajectory_XYZ(j,i*3-1) = cosd(phi) * sind(lambda); %y position
        trajectory_XYZ(j,i*3) = sind(phi);
    end
    
    
end

