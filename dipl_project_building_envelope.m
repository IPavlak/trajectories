x = [10 15 20 30 10];
y = [10 18 10 30 20];
x = [10 50 50 40 40 20 20 10 ];
y = [10 10 30 30 20 20 30 30 ];
res = 1;
res_building_points = 0.5;
dist = 2;

[mx my] = weight_center(x, y);

plot(x, y, '.', 'MarkerSize', 20);
xlim([0, 30]); 
ylim([0, 30]);

[x, y] = interpolate(x, y, res_building_points);

plot(x, y, 'k.', 'MarkerSize', 20);
xlim([0, 40]); 
ylim([0, 40]);
hold on;
axis equal;

[traj_x, traj_y] = intersection(x(1), y(1), x(end), y(end), dist, dist, mx, my, 'd');

s = 0;
start_inter_x = traj_x;
start_inter_y = traj_y;
i = 1;
minAngleIndex = 1;
while minAngleIndex >= i
    i = minAngleIndex;
    minAngle = 3*pi;
    
    cnt = 0; j = i; check_pts = length(x)/4;
    end_inter_x = 0; end_inter_y = 0;
    while cnt < check_pts
%         j = j+1;
        j = mod(j, length(x))+1;
        cnt = cnt+1;
        
        if sqrt( (x(i)-x(j))^2 + (y(i)-y(j))^2 ) > 2*dist
            continue;
        end
        
        [inter_x, inter_y] = intersection_ccw(x(i), y(i), x(j), y(j), dist, dist, start_inter_x, start_inter_y);
        angle = getAngle_ccw(start_inter_x, start_inter_y, inter_x, inter_y, x(i), y(i));
        if(angle < minAngle)
            minAngle = angle;
            minAngleIndex = j;
            end_inter_x = inter_x; end_inter_y = inter_y;
        end
    end
    start_inter_x = end_inter_x;
    start_inter_y = end_inter_y;
end
i = minAngleIndex;
traj_x = end_inter_x;
traj_y = end_inter_y;

while minAngleIndex >= i
    i = minAngleIndex;
% for i=1:13%length(x)
    minAngle = 3*pi;
    
    cnt = 0; j = i; check_pts = length(x)/4;
    end_inter_x = 0; end_inter_y = 0;
    while cnt < check_pts
%         j = j+1;
        j = mod(j, length(x))+1;
        cnt = cnt+1;
        
        if sqrt( (x(i)-x(j))^2 + (y(i)-y(j))^2 ) > 2*dist
            continue;
        end
        
        [inter_x, inter_y] = intersection_ccw(x(i), y(i), x(j), y(j), dist, dist, start_inter_x, start_inter_y);
        angle = getAngle_ccw(start_inter_x, start_inter_y, inter_x, inter_y, x(i), y(i));
        if(angle < minAngle)
            minAngle = angle;
            minAngleIndex = j;
            end_inter_x = inter_x; end_inter_y = inter_y;
        end
    end
%     plot(end_inter_x, end_inter_y, 'r.', 'MarkerSize', 20);
    
    [new_traj_points_x, new_traj_points_y, s] = ...
        make_traj_pts(start_inter_x, start_inter_y, end_inter_x, end_inter_y, x(i), y(i), s, res, dist);

    traj_x = [traj_x, new_traj_points_x];
    traj_y = [traj_y, new_traj_points_y];
    
    start_inter_x = end_inter_x;
    start_inter_y = end_inter_y;
end

plot(traj_x, traj_y, 'b.', 'MarkerSize', 15);

% for i=1:length(x)
%     circle(x(i), y(i), dist)
% end



function [new_xs, new_ys, new_s] = make_traj_pts(x1, y1, x2, y2, sx, sy, s, res, dist)
new_xs = []; new_ys = [];
% circle(sx, sy, dist);
% plot([x1 x2], [y1 y2], 'r.', 'MarkerSize', 15);
angle = getAngle_ccw(x1, y1, x2, y2, sx, sy);
l = dist*angle;
while l > arc_over_secant(res-s, dist)
%     [x1, y1] = intersection(sx, sy, x1, y1, dist, res-s, x2, y2, 'b');
    [x1, y1] = intersection_ccw(sx, sy, x1, y1, dist, res-s, x1, y1);
    new_xs = [new_xs, x1]; new_ys = [new_ys, y1];
    
%     angle = getAngle(x1, y1, x2, y2, dist);
    angle = getAngle_ccw(x1, y1, x2, y2, sx, sy);
    l = dist*angle;
    s = 0;
end
new_s = s+l;
end



function [interp_x, interp_y] = interpolate(x, y, res)
interp_x = [];
interp_y = [];
for i=1:length(x)-1
    n = sqrt( (x(i+1)-x(i))^2 + (y(i+1)-y(i))^2 ) / res;
    interp_x = [interp_x, linspace(x(i), x(i+1), n)];
    interp_y = [interp_y, linspace(y(i), y(i+1), n)];
end

    n = sqrt( (x(1)-x(end))^2 + (y(1)-y(end))^2 ) / res;
    interp_x = [interp_x, linspace(x(end), x(1), n)];
    interp_y = [interp_y, linspace(y(end), y(1), n)];

    
del_index = [];
for i=2:length(interp_x)
    if interp_x(i) == interp_x(i-1) && interp_y(i) == interp_y(i-1)
        del_index = [del_index, i];
    end
end
if interp_x(end) == interp_x(1) && interp_y(end) == interp_y(1)
    del_index = [del_index, length(interp_x)];
end
interp_x(del_index) = [];
interp_y(del_index) = [];

end



function [x, y] = intersection(sx1, sy1, sx2, sy2, r1, r2, mx ,my, cond)
R = sqrt( (sx2-sx1)^2 + (sy2-sy1)^2 );

x1 = 1/2*(sx1+sx2) + (r1^2-r2^2)/(2*R^2)*(sx2-sx1) + 1/2*sqrt(2*(r1^2+r2^2)/R^2 - (r1^2-r2^2)^2/R^4 - 1)*(sy2-sy1);
x2 = 1/2*(sx1+sx2) + (r1^2-r2^2)/(2*R^2)*(sx2-sx1) - 1/2*sqrt(2*(r1^2+r2^2)/R^2 - (r1^2-r2^2)^2/R^4 - 1)*(sy2-sy1);

y1 = 1/2*(sy1+sy2) + (r1^2-r2^2)/(2*R^2)*(sy2-sy1) + 1/2*sqrt(2*(r1^2+r2^2)/R^2 - (r1^2-r2^2)^2/R^4 - 1)*(sx1-sx2);
y2 = 1/2*(sy1+sy2) + (r1^2-r2^2)/(2*R^2)*(sy2-sy1) - 1/2*sqrt(2*(r1^2+r2^2)/R^2 - (r1^2-r2^2)^2/R^4 - 1)*(sx1-sx2);

% circle(sx1,sy1,r1); hold on;
% circle(sx2,sy2,r2);
% plot([x1 x2], [y1 y2], 'b.', 'MarkerSize', 10);
% xlim([4 18]); ylim([4 16]);
% getAngle(x1, y1, x2, y2, sx1, sx2, dist)

d1 = sqrt( (x1-mx)^2 + (y1-my)^2 );
d2 = sqrt( (x2-mx)^2 + (y2-my)^2 );

if(d1 > d2)
    if cond == 'd'
        x = x1; y = y1;
    else
        x = x2; y = y2;
    end
else
    if cond == 'b'
        x = x1; y = y1;
    else
        x = x2; y = y2;
    end
end
% plot([x], [y], 'b.', 'MarkerSize', 10);
end


function [x_inter, y_inter] = intersection_ccw(sx1, sy1, sx2, sy2, r1, r2, x, y)
R = sqrt( (sx2-sx1)^2 + (sy2-sy1)^2 );

x1 = 1/2*(sx1+sx2) + (r1^2-r2^2)/(2*R^2)*(sx2-sx1) + 1/2*sqrt(2*(r1^2+r2^2)/R^2 - (r1^2-r2^2)^2/R^4 - 1)*(sy2-sy1);
x2 = 1/2*(sx1+sx2) + (r1^2-r2^2)/(2*R^2)*(sx2-sx1) - 1/2*sqrt(2*(r1^2+r2^2)/R^2 - (r1^2-r2^2)^2/R^4 - 1)*(sy2-sy1);

y1 = 1/2*(sy1+sy2) + (r1^2-r2^2)/(2*R^2)*(sy2-sy1) + 1/2*sqrt(2*(r1^2+r2^2)/R^2 - (r1^2-r2^2)^2/R^4 - 1)*(sx1-sx2);
y2 = 1/2*(sy1+sy2) + (r1^2-r2^2)/(2*R^2)*(sy2-sy1) - 1/2*sqrt(2*(r1^2+r2^2)/R^2 - (r1^2-r2^2)^2/R^4 - 1)*(sx1-sx2);

angle1 = getAngle_ccw(x, y, x1, y1, sx1, sy1);
angle2 = getAngle_ccw(x, y, x2, y2, sx1, sy1);

% circle(sx1,sy1,r1); hold on;
% circle(sx2,sy2,r2);
% plot([x1 x2], [y1 y2], 'b.', 'MarkerSize', 10);
% xlim([4 18]); ylim([4 16]);
% getAngle(x1, y1, x2, y2, sx1, sx2, dist)

if angle1 < angle2
    x_inter = x1; y_inter = y1;
else
    x_inter = x2; y_inter = y2;
end
% plot([x_inter], [y_inter], 'r.', 'MarkerSize', 10);
end


function alpha = getAngle(x1, y1, x2, y2, R)
t = sqrt( (x1-x2)^2 + (y1-y2)^2 );  %tetive
cos_alpha = (2*R^2 - t^2)/(2*R^2);
alpha = acos(cos_alpha);
end


function alpha = getAngle_ccw(x1, y1, x2, y2, sx, sy)
x1 = x1 - sx; y1 = y1 - sy;
x2 = x2 - sx; y2 = y2 - sy;

alpha1 = atan2(y1, x1);
alpha2 = atan2(y2, x2);

alpha = alpha2 - alpha1;
if alpha < 0
    alpha = 2*pi + alpha;
end
end


function x = arc_over_secant(d, r)
alpha = acos( (2*r^2 - d^2) / (2*r^2) );
x = r*alpha;
end



function [mx, my] = weight_center(x, y)
sumx = sum(x);
sumy = sum(y);
mx = sumx/length(x);
my = sumy/length(y);
end


function circle(x,y,r)
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
plot(xunit, yunit, 'g.');
end