clear
close all


data = csvread('ThirtyRobotsAsteroid_ellipsoid.csv');
objectType = 'asteroid_ellipsoid';
% objectType can be 'flat_earth'  or 'rectangular_prism'
r_circle=5;
n = 30; %number of robots
dt = 0.1;
[length, width] = size(data);
timesteps = length/n; %number of timesteps
for i=1:1:timesteps
    %trajectory_angles is an n by 2*timesteps matrix containing time history of
    %all robots
    current_angles = data(i*n-(n-1):i*n,:);
    trajectory_angles(1:n,i*2-1:i*2) = current_angles;
    for j=1:1:n
    %get trajectory_XYZ, a n by 3*timestep matrix containing time history
    %of all robots in XYZ coordinates
        phi = current_angles(j,1);
        lambda = current_angles(j,2);
        trajectory_XYZ(j,i*3-2) = r_circle*cosd(phi) * cosd(lambda); %x position
        trajectory_XYZ(j,i*3-1) = r_circle*cosd(phi) * sind(lambda); %y position
        trajectory_XYZ(j,i*3) = r_circle*sind(phi);
    end
    
    
end

%Animation
f = figure('Color', [1 1 1 1], 'Position', [403 176 698 490])
if (strcmp(objectType, 'rectangular_prism'))
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
elseif (strcmp(objectType,'flat_earth'))
    %% make flat plate
    xc=0; yc=0; zc=0;    % coordinated of the center
    L=0.6;                 % width and length size (length of an edge)
    H=L/20;
    alpha=0.3;           % transparency (max=1=opaque)
    X = [0 0 0 0 0 1; 1 0 1 1 1 1; 1 0 1 1 1 1; 0 0 0 0 0 1];
    Y = [0 0 0 0 1 0; 0 1 0 0 1 1; 0 1 1 1 1 1; 0 0 1 1 1 0];
    Z = [0 0 1 0 0 0; 0 0 1 0 0 0; 1 1 1 0 1 1; 1 1 1 0 1 1];
    C='blue';                  % unicolor
    X = L*(X-0.5) + xc;
    Y = L*(Y-0.5) + yc;
    Z = H*(Z-0.5) + zc; 
    
elseif (strcmp(objectType,'asteroid_ellipsoid'))
    %% make ellipsoid asteroid
    theta = gallery('uniformdata',[100,1],0)*2*pi;
    phi = gallery('uniformdata',[100,1],1)*pi;
    x = 2.*cos(theta).*sin(phi);
    y = sin(theta).*sin(phi);
    z = cos(phi);
    DT = delaunayTriangulation(x,y,z);
    [T,Xb] = freeBoundary(DT);
    TR = triangulation(T,Xb);
    P = incenter(TR);
    F = faceNormal(TR);  
end


%% make circle equator and prime meridians
teta=-pi:0.01:pi;
psiRange=-pi:pi/12:pi;
x_equator=r_circle*cos(teta);
y_equator=r_circle*sin(teta);
x_meridian=r_circle*cos(teta);
z_meridian=r_circle*sin(teta);
y_extra=r_circle*cos(teta);
z_extra=r_circle*sin(teta);
% x_others=[]
for j=1:size(psiRange,2)
    psi=psiRange(j);
    x_others(j,:) = r_circle*cos(teta)*sin(psi);
    y_others(j,:) = r_circle*cos(teta)*cos(psi);
    z_others(j,:) = r_circle*sin(teta);
end
%%
set(f, 'doublebuffer', 'on');
max_t = timesteps; %Find end time
t = 0;  %Set movie time to 0
i = 1;  %Set index of array to start
pause(1);

v = VideoWriter('Space_Force_Two.avi');
v.FrameRate = 1/0.1;

while i<=timesteps+50
    %Every dt seconds, show position of robots
    %frame
    if i<=timesteps
        for j=1:1:n
            % add robots of this timeframe
            plot3(trajectory_XYZ(j,i*3-2),trajectory_XYZ(j,i*3-1),trajectory_XYZ(j,i*3), 'o','MarkerFaceColor',[1 .6 .6])  
            if j ==1
                hold on
            end
        end
    elseif i<=timesteps+15
        for j=1:1:n
            % add robots of this timeframe
            plot3(trajectory_XYZ(j,timesteps*3-2),trajectory_XYZ(j,timesteps*3-1),trajectory_XYZ(j,timesteps*3), 'o','MarkerFaceColor',[1 .6 .6])
            hold on
            plot3([trajectory_XYZ(j,timesteps*3-2) 0],[trajectory_XYZ(j,timesteps*3-1) 0],[trajectory_XYZ(j,timesteps*3) 0], 'r','MarkerFaceColor',[1 .6 .6])  
            if j ==1
                hold on
            end
        end
    else
        for j=1:1:n
            % add robots of this timeframe
            plot3(trajectory_XYZ(j,timesteps*3-2),trajectory_XYZ(j,timesteps*3-1),trajectory_XYZ(j,timesteps*3), 'o','MarkerFaceColor',[1 .6 .6]) 
            if j ==1
                hold on
            end
        end
    end
    
    %fill3(X,Y,Z,C,'FaceAlpha',alpha) %add cube
    if i<=timesteps+15
        trisurf(T,Xb(:,1),Xb(:,2),Xb(:,3), ...
        'FaceColor','cyan','FaceAlpha',0.8);
    elseif i<=timesteps+19
        for j=1:1:n
            % add robots of this timeframe
            img = imread('explosion.png'); %add cube
            image('CData',img,'XData',[-2 2],'YData',[-2 2])
            
            if j ==1
                hold on
            end
        end
    else
        test = 0;
    end
    
    plot3(x_equator,y_equator,zeros(1,numel(x_equator)),':r', 'LineWidth', 0.5)
   
    for j=1:size(psiRange,2)
        plot3(x_others(j,:),y_others(j,:),z_others(j,:),':b','LineWidth',0.5)
    end
    hold off
    axis equal
    %axis([-1.2, 1.2, -1.2, 1.2, -1.2, 1.2])
    ax1 = gca;
    ax1.XAxis.Visible = 'off';
    ax1.YAxis.Visible = 'off';
    ax1.ZAxis.Visible = 'off';
    
    %view(3)
    view(i*(120/timesteps),20)
    xlabel('x axis')
    ylabel('y axis')
    if i<=timesteps
        title(strcat('Time = ', num2str(t,2),' seconds'))
    else
        title('Converged')
    end
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

