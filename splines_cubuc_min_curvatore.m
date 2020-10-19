x0 = 200; y0 = 200; theta0 = -130*pi/180;
xf = 300; yf = 100;

epsilon = 1e-5;     % when are 2 values equal
N = 1000;           % number of points between (x0,y0) and (xf,yf)

%% minimum curvature (cubic splines)
ka = xf - x0;
kb = yf - y0;
tg0 = tan(theta0)

if abs(tg0) > 10
    ctg0 = cot(theta0);
    A = [1, 3, 0, 0, 0;
         0, 0, 0, 1, 3;
         ctg0, 0, 0, 1, 0;
         1, 1, ctg0, 0, 0;
         0, 0, 1, 1, 1];
else
    A = [0, 1, 3, 0, 0;
         0, 0, 0, 1, 3;
         0, 1, 0, tg0, 0;
         1, 1, 1, 0, 0;
         tg0, 0, 0, 1, 1 ];
end
 
inv_A = inv(A);
spline_coeff = inv_A * [0, 0, 0, ka, kb]'

if abs(tg0) > 10
    a0 = x0; a2 = spline_coeff(1); a3 = spline_coeff(2); 
    b0 = y0; b1 = spline_coeff(3); b2 = spline_coeff(4); b3 = spline_coeff(5);
    a1 = b1*ctg0;
else
    a0 = x0; a1 = spline_coeff(1); a2 = spline_coeff(2); a3 = spline_coeff(3);
    b0 = y0; b1 = a1*tg0; b2 = spline_coeff(4); b3 = spline_coeff(5);
end

% wrong angle
if abs(atan2(b1, a1) - theta0) > epsilon
    disp('tu sam')
    A = [1, 3, 0, 0;
         0, 0, 1, 3;
         1, 1, 0, 0;
         0, 0, 1, 1];
     inv_A = inv(A);
     spline_coeff = inv_A * [0, 0, ka+a1, kb+b1]'


    if abs(tg0) > 10
        a0 = x0; a2 = spline_coeff(1); a3 = spline_coeff(2); 
        b0 = y0; b1 = -b1; b2 = spline_coeff(3); b3 = spline_coeff(4);
        a1 = b1*ctg0;
    else
        a0 = x0; a1 = -a1; a2 = spline_coeff(1); a3 = spline_coeff(2);
        b0 = y0; b1 = a1*tg0; b2 = spline_coeff(3); b3 = spline_coeff(4);
    end
end

% a0 = x0; a1 = 1000*cos(theta0);
% b0 = y0; b1 = 1000*sin(theta0);
% ka = xf - a0 - a1; kb = yf - b0 - b1;
% a3 = -ka/2; b3 = -kb/2;
% a2 = -3*a3; b2 = -3*b3;
% a = [a3 a2 a1 a0]
% b = [b3 b2 b1 b0]

u = linspace(0, 1, N);
x = a3*u.^3 + a2*u.^2 + a1*u + a0;
y = b3*u.^3 + b2*u.^2 + b1*u + b0;

plot(x,y)

ddx = 6*a3*u + 2*a2;
ddy = 6*b3*u + 2*b2;
% plot(u, ddy)

%% minimum curvature (quintic splines)
a0 = x0; b0 = y0;
a1 = 140*cos(theta0); b1 = 140*sin(theta0); % K = 140
ka = xf - a0 - a1; kb = yf - b0 - b1;

A = [ 2 6 12 20 0 0 0 0 1 0;
      6 18 24 60 0 0 0 0 1 0;
      12 24 48 120 0 0 0 0 1 0;
      20 60 120 200 0 0 0 0 1 0;
      0 0 0 0 2 6 12 20 0 1;
      0 0 0 0 6 18 24 60 0 1;
      0 0 0 0 12 24 48 120 0 1;
      0 0 0 0 20 60 120 200 0 1;
      1 1 1 1 0 0 0 0 0 0;
      0 0 0 0 1 1 1 1 0 0           ];
  
inv_A = inv(A);
% spline_coeff = inv_A * [0 0 0 0 0 0 0 0 ka kb]';
% a = [a0 a1 spline_coeff(1) spline_coeff(2) spline_coeff(3) spline_coeff(4)]
% b = [b0 b1 spline_coeff(5) spline_coeff(6) spline_coeff(7) spline_coeff(8)]
[a, b, k] = bestK(a0, b0, xf, yf, theta0, inv_A);
u = linspace(0, 1, N);
x = a(6)*u.^5 + a(5)*u.^4 + a(4)*u.^3 + a(3)*u.^2 + a(2)*u + a(1);
y = b(6)*u.^5 + b(5)*u.^4 + b(4)*u.^3 + b(3)*u.^2 + b(2)*u + b(1);

ddx = 20*a(6)*u.^3 + 12*a(5)*u.^2 + 6*a(4)*u + 2*a(3);
ddy = 20*b(6)*u.^3 + 12*b(5)*u.^2 + 6*b(4)*u + 2*b(3);
trapz(u, ddx.*ddx);
figure(1);
plot(x,y)

%% minimum curvature -> starting curvature = 0
a0 = x0; b0 = y0;
a1 = 100*cos(theta0); b1 = 100*sin(theta0);
ka = xf - a0 - a1; kb = yf - b0 - b1;

A = [ 18 24 60 0 0 0 1 0 6 0;
      24 48 120 0 0 0 1 0 12 0;
      60 120 200 0 0 0 1 0 20 0;
      0 0 0 18 24 60 0 1 0 6;
      0 0 0 24 48 120 0 1 0 12;
      0 0 0 60 120 200 0 1 0 20;
      1 1 1 0 0 0 0 0 0 0;
      0 0 0 1 1 1 0 0 0 0;
      6 12 20 0 0 0 0 0 0 0;
      0 0 0 6 12 20 0 0 0 0];

inv_A = inv(A);
spline_coeff = inv_A * [0 0 0 0 0 0 ka kb 0 0]';
a = [a0 a1  0 spline_coeff(1) spline_coeff(2) spline_coeff(3)]
b = [b0 b1  0 spline_coeff(4) spline_coeff(5) spline_coeff(6)]

u = linspace(0, 1, N);
x = a(6)*u.^5 + a(5)*u.^4 + a(4)*u.^3 + a(3)*u.^2 + a(2)*u + a(1);
y = b(6)*u.^5 + b(5)*u.^4 + b(4)*u.^3 + b(3)*u.^2 + b(2)*u + b(1);

ddx = 20*a(6)*u.^3 + 12*a(5)*u.^2 + 6*a(4)*u + 2*a(3);
ddy = 20*b(6)*u.^3 + 12*b(5)*u.^2 + 6*b(4)*u + 2*b(3);
trapz(u, ddx.*ddx)
plot(u,ddy)

%% constant curvature
tg0 = tan(theta0);

% special cases
if abs(atan2(yf-y0, xf-x0) - theta0) < epsilon
    disp('colinear points')
    x = linspace(x0, xf, N);
    y = linspace(y0, yf, N);
    return;
elseif abs(atan2(yf-y0, xf-x0) + theta0) < epsilon
    disp('colienar points - wrong direction (not allowed)')
    return;

elseif abs(tg0) > 1e9
    disp("vertical gradient")
    k1 = (x0^2 + y0^2 - xf^2 - yf^2) / (2*(x0-xf));
    k2 = (y0-yf) / (x0-xf);
    cy = y0
    cx = k1 - cy*k2
elseif x0 == xf
    disp("x0 == xf")
    cy = (y0+yf) / 2
    cx = x0 - tg0*(cy-y0)
    if abs( (cx-x0) / (cy-y0) + tg0) > epsilon % maybe unneccesary
        disp('wrong angle')
    end
else
    k1 = (x0^2 + y0^2 - xf^2 - yf^2) / (2*(x0-xf));
    k2 = (y0-yf) / (x0-xf);
    cy = (-y0*tg0 + k1 - x0) / (k2 - tg0)
    cx = k1 - cy*k2
    if abs((cx-x0) / (cy-y0) + tg0) > epsilon
        disp(['wrong angle: ', num2str((cx-x0)/(cy-y0)), num2str(tg0)])
        cy = (y0*tg0 + k1 - x0) / (k2 + tg0)
        cx = k1 - cy*k2
    end
end

R = sqrt((x0-cx)^2 + (y0-cy)^2)


% circle direction
alpha = atan2((y0-cy), (x0-cx))
betta = atan2((yf-cy), (xf-cx))

D = [x0-cx, y0-cy; cos(theta0), sin(theta0)];
deter = det(D);
if deter > 0
    if alpha == pi
        alpha = -pi;
    end
    if betta < alpha
        betta = betta + 2*pi;
    end
    thetaf = theta0 + (betta-alpha)
    u = linspace(alpha, betta, N);
else
    if alpha == -pi
        alpha = pi;
    end
    if betta > alpha
        betta = betta - 2*pi;
    end
    thetaf = theta0 - (alpha-betta)
    u = linspace(alpha, betta, N);
end

x = cx + R*cos(u);
y = cy + R*sin(u);

plot(x, y);


%%
function [a, b, K] = bestK(a0, b0, xf, yf, theta0, inv_A)
N = 100;
u = linspace(0, 1, N);
delta_theta = zeros(1, N-1);
ks = [[0.1:0.1:1] [1:200] [200:10:1000] [1000:100:10000]];
n = size(ks, 2);
curveInts = zeros(1, n);
for j = 1:n
    k = ks(j);
    a1 = k*cos(theta0);
    b1 = k*sin(theta0);
    ka = xf - a0 - a1; kb = yf - b0 - b1;
    b = [0 0 0 0 0 0 0 0 ka kb]';
    spline_coeff = inv_A * b;
    a = [a0 a1 spline_coeff(1) spline_coeff(2) spline_coeff(3) spline_coeff(4)];
    b = [b0 b1 spline_coeff(5) spline_coeff(6) spline_coeff(7) spline_coeff(8)];
    x = a(6)*u.^5 + a(5)*u.^4 + a(4)*u.^3 + a(3)*u.^2 + a(2)*u + a(1);
    y = b(6)*u.^5 + b(5)*u.^4 + b(4)*u.^3 + b(3)*u.^2 + b(2)*u + b(1);
    
    delta_theta(1) = angdiff( atan2(y(2)-y(1), x(2)-x(1)), theta0)^2;
    if delta_theta(1) > 0.8
        delta_theta(1) = Inf;
    end
    for i = 2:N-1
        delta_theta(i) = angdiff( atan2(y(i)-y(i-1), x(i)-x(i-1)), atan2(y(i+1)-y(i), x(i+1)-x(i)) )^2;
        if delta_theta(i) > 0.8
            delta_theta(i) = Inf;
        end
    end
    curveInts(j) = trapz(delta_theta);
    k
end
figure(2);
plot(curveInts);
[M, I] = min(curveInts);
M
K = ks(I)
a1 = K*cos(theta0);
b1 = K*sin(theta0);
ka = xf - a0 - a1; kb = yf - b0 - b1;
b = [0 0 0 0 0 0 0 0 ka kb]';
spline_coeff = inv_A * b;
a = [a0 a1 spline_coeff(1) spline_coeff(2) spline_coeff(3) spline_coeff(4)];
b = [b0 b1 spline_coeff(5) spline_coeff(6) spline_coeff(7) spline_coeff(8)];

% x = a(6)*u.^5 + a(5)*u.^4 + a(4)*u.^3 + a(3)*u.^2 + a(2)*u + a(1);
% y = b(6)*u.^5 + b(5)*u.^4 + b(4)*u.^3 + b(3)*u.^2 + b(2)*u + b(1);
% delta_theta(1) = angdiff( atan2(y(2)-y(1), x(2)-x(1)), theta0)^2;
% for i = 2:N-1
%     delta_theta(i) = angdiff( atan2(y(i)-y(i-1), x(i)-x(i-1)), atan2(y(i+1)-y(i), x(i+1)-x(i)) )^2;
% end
% plot(u(1:N-1), delta_theta);
end