%Ritvik Verma
%Project 2
%Phase 1
%Modelling a home run hit with air resistance
%Comparing the analytic and the numeric solutions without drag

clear
clf

v0mph = 112;   % exit velocity in mph 
phi0deg = 32;    % launch angle in degrees 

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

m = 0.145;      %mass of the baseball in kg

dt = (tmax-tmin)/N;
y = zeros(1, N+1);   % initialize y(t)

y(1) = y0;
x(1) = x0;
vy = v0y;       %setting initial velocities and positions for x and y
vx = v0x;

for n = 1:N   % stop at N
    Fnet_x = 0;     
    Fnet_y = -m*g; 
    ax = Fnet_x/m;
    ay = Fnet_y/m;   
    y(n+1) = y(n) + vy*dt + (1/2)*ay*dt^2;
    vy = vy + ay*dt;
    x(n+1) = x(n) + vx*dt + (1/2)*ax*dt^2;
    vx = vx + ax*dt;
end

x_ft = x*m2ft;
y_ft = y*m2ft;

checkSum_y = sum(abs(y_ft-yt)) 
checkSum_x = sum(abs(x_ft-xt)) 

%---------------plotting the trajectory---------------

plot(xt, yt, x_ft, y_ft,'LineWidth', 2) 
grid on
ax = gca; ax.FontSize = 16; 
ax.GridAlpha = 0.4;

title({'ECE 202, Project 2, Phase 1: Trajectory of a baseball with no drag', ...
    'Analytic vs. Numeric solution'}, 'FontSize', 22)
xlabel('x (ft)', 'FontSize', 18)   
ylabel('y (ft)', 'FontSize', 18)
legend({'analytic (behind numeric)', 'numeric'}, ...
    'FontSize', 18)
