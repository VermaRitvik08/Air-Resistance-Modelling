%Ritvik Verma
%Project 2
%Phase 3
%Modelling a home run hit with air resistance
%Exploring data and analyzing it in Excel

clear
clf

v0mph = 112;   % exit velocity in mph 
phi0deg = 32;    % launch angle in degrees 
C = input('Enter value for C: ');
m = 0.145;      %mass of the baseball in kg
A = 0.0042;     %cross-section area of baseball units in m^2
p = 1.225;      %density of air units in kg/m^3


x0 = 0;         %start coordinates of ball 
y0 = 0;
g = 10;     % gravitational constant in N/kg

mph2mps = 5280 * 12 * 2.54 / 100 / 3600;   % mph to m/s conversion
deg2rad = pi()/180;   % degrees to radians

v0 = v0mph * mph2mps;  
phi0 = phi0deg * deg2rad;

v0x = v0*cos(phi0);   % x-component of v0
v0y = v0*sin(phi0);   % y-component of v0

% ----- compute some useful characteristics of trajectory -----

tH = v0y/g;    % time to reach max. height
tLand = 2*tH;   % time to land (time of flight)

H = tH * v0y/2;   % max. height
R = v0x * tLand;   % range

m2ft = 3.28084;     %conversion constant from meter to ft

% ----- set up a time array, compute x(t), y(t) analytically -----

tmin = 0; 
tmax = tLand; 
N = 2000;   % intervals

t = linspace(tmin, tmax, N+1);   % time array, connects x(t) with y(t)

xt = (x0 + v0x*t)*m2ft;   
yt = (y0 + v0y*t - (1/2)*g*t.^2)*m2ft;   

% ----- numeric solution -----

dt = (tmax-tmin)/N;
y = zeros(1, N+1);   % initialize y(t)

y(1) = y0;
x(1) = x0;
vy = v0y;       %setting initial velocities and positions for x and y
vx = v0x;

D = 0.5*C*p*A; %positive constant in drag force

for n = 1:N   % stop at N
    v = sqrt(vx^2 + vy^2);
    Fnet_x = 0 - D*vx*v;     
    Fnet_y = -m*g - D*vy*v; 
    ax = Fnet_x/m;
    ay = Fnet_y/m;   
    y(n+1) = y(n) + vy*dt + (1/2)*ay*dt^2;
    vy = vy + ay*dt;
    x(n+1) = x(n) + vx*dt + (1/2)*ax*dt^2;
    vx = vx + ax*dt;
end

x_ft = x*m2ft;
y_ft = y*m2ft;

if C == 0
    checkSum_y = sum(abs(y_ft-yt)) 
    checkSum_x = sum(abs(x_ft-xt))
end
%---------------plotting the trajectory---------------

plot(xt, yt, x_ft, y_ft,'LineWidth', 2) 
grid on
set(gca,'XMinorGrid','on');
set(gca,'YMinorGrid','on');
ax = gca; ax.FontSize = 16; 
ax.GridAlpha = 0.4;
ax.MinorGridAlpha = 0.5;

ylim([0 150])
title({'ECE 202, Project 2, Phase 3: Trajectory of a baseball', ...
    'drag vs without drag'}, 'FontSize', 22)
xlabel('x (ft)', 'FontSize', 18)   
ylabel('y (ft)', 'FontSize', 18)
str1 = sprintf('with drag and C = %g',C);
legend({'without drag', str1}, ...
    'FontSize', 18)

export = [t; x_ft; y_ft].';
labels = ["Time (s)", "x (ft)", "Height (ft)"];

export = [labels; export];
writematrix(export, 'p3_2.csv')


