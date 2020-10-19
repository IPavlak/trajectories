axis equal

P = rand(30, 3);
% [x, y, z] = spherePoints(0, 0, 0, 1);
% P = [x, y, z];

plot3(P(:,1),P(:,2),P(:,3),'.','MarkerSize',10)
grid on

k = boundary(P);
trisurf(k,P(:,1),P(:,2),P(:,3),'Facecolor','red','FaceAlpha',0.1)
xlim([-3, 3])
ylim([-3, 3])
zlim([-3, 3])
hold on

for i=1:size(P, 1)
    [x,y,z] = spherePoints(P(i,1), P(i,2), P(i,3), 1);
    X = [P(:,1); x];
    Y = [P(:,2); y];
    Z = [P(:,3); z];
    P = [X, Y, Z];
end

% plot3(P(:,1),P(:,2),P(:,3),'.','MarkerSize',10)
k = boundary(P);
trisurf(k,P(:,1),P(:,2),P(:,3),'Facecolor','blue','FaceAlpha',0.1)


%%
% spherePoints(0, 0, 0, 1)

function [x, y, z] = spherePoints(cx, cy, cz, r)
    theta=linspace(0,2*pi,10);
    phi=linspace(0,pi,10);
    [theta,phi]=meshgrid(theta,phi);
    rho=1;
    x=rho*sin(phi).*cos(theta);
    y=rho*sin(phi).*sin(theta);
    z=rho*cos(phi);
%     mesh(x,y,z)
%     plot3(2*x(:)+1, 2*y(:)+1, 2*z(:), '.')
    x = r*x(:) + cx;
    y = r*y(:) + cy;
    z = r*z(:) + cz;
end