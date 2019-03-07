clear
close all

data = xlsread('Positions.csv');
n = 6; %number of robots
dt = 0.1;
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

%Animation
f = figure('Color', [1 1 1 1], 'Position', [403 176 698 490])
%% make cube
xc=0; yc=0; zc=0;    % coordinated of the center
L=0.6;                 % cube size (length of an edge)
alpha=0.3;           % transparency (max=1=opaque)
X = [0 0 0 0 0 1; 1 0 1 1 1 1; 1 0 1 1 1 1; 0 0 0 0 0 1];
Y = [0 0 0 0 1 0; 0 1 0 0 1 1; 0 1 1 1 1 1; 0 0 1 1 1 0];
Z = [0 0 1 0 0 0; 0 0 1 0 0 0; 1 1 1 0 1 1; 1 1 1 0 1 1];
C='blue';                  % unicolor
X = L*(X-0.5) + xc;
Y = L*(Y-0.5) + yc;
Z = L*(Z-0.5) + zc; 

%% make circle equator and prime meridians
r_circle=1;
teta=-pi:0.01:pi;
x_equator=r_circle*cos(teta);
y_equator=r_circle*sin(teta);
x_meridian=r_circle*cos(teta);
z_meridian=r_circle*sin(teta);


%%
set(f, 'doublebuffer', 'on');
max_t = timesteps; %Find end time
t = 0;  %Set movie time to 0
i = 1;  %Set index of array to start
pause(1);

v = VideoWriter('Space_Force_One.avi');
v.FrameRate = 1/0.1;

while i<=timesteps
    %Every dt seconds, show position of robots
    %frame
    
    for j=1:1:n
        % add robots of this timeframe
        plot3(trajectory_XYZ(j,i*3-2),trajectory_XYZ(j,i*3-1),trajectory_XYZ(j,i*3), 'o','MarkerFaceColor',[1 .6 .6])  
        if j ==1
            hold on
        end
    end
    fill3(X,Y,Z,C,'FaceAlpha',alpha) %add cube
    plot3(x_equator,y_equator,zeros(1,numel(x_equator)))
    plot3(x_meridian,zeros(1,numel(x_meridian)),z_meridian)
    hold off
    axis equal
    axis([-1.2, 1.2, -1.2, 1.2, -1.2, 1.2])
    view(280,10)
    xlabel('x axis')
    ylabel('y axis')
    title(strcat('Time = ', num2str(t,2),' seconds'))
    grid minor
    drawnow;
    M(i) = getframe(1);
    pause(0.1)
    t = t+dt;
    i=i+1;
end

open(v);
writeVideo(v,M);
close(v)

