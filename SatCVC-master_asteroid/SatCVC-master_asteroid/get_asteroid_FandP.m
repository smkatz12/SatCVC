
clc;
clear;

theta = gallery('uniformdata',[100,1],0)*2*pi;
phi = gallery('uniformdata',[100,1],1)*pi;
%x = cos(theta).*sin(phi);
%y = sin(theta).*sin(phi);
%z = cos(phi);
x = 2.*cos(theta).*sin(phi);
y = sin(theta).*sin(phi);
z = cos(phi);




DT = delaunayTriangulation(x,y,z);
[T,Xb] = freeBoundary(DT);
% for i = 1:2:length(Xb)
%     if Xb(i,1) < 0
%         Xb(i,1) = Xb(i,1)-0.1;
%     else
%         Xb(i,1) = Xb(i,1)+0.1;
%     end
%     if Xb(i,2) < 0
%         Xb(i,2) = Xb(i,2)-0.1;
%     else
%         Xb(i,2) = Xb(i,2)+0.1;
%     end
%     if Xb(i,3) < 0
%         Xb(i,3) = Xb(i,3)-0.1;
%     else
%         Xb(i,3) = Xb(i,3)+0.1;
%     end
% end

TR = triangulation(T,Xb);
P = incenter(TR);
F = faceNormal(TR);  
trisurf(T,Xb(:,1),Xb(:,2),Xb(:,3), ...
     'FaceColor','cyan','FaceAlpha',0.8);
axis equal
hold on  
quiver3(P(:,1),P(:,2),P(:,3), ...
     F(:,1),F(:,2),F(:,3),0.5,'color','r');
 
csvwrite('F_asteroid_ellipsoid.csv',F)
csvwrite('P_asteroid_ellipsoid.csv',P)

