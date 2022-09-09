%Ritvik Verma
%Project 2
%Phase 5
%Modelling a home run hit with air resistance
%Exploring results 

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
    
    if y(n)/y(n+1) <= 0
        time_of_flight_s = t(n); 
        xf = x(n);  %distance where the ratio is negative i.e when ball hit
                    %the ground
        vf = v;     %velocity when the ball approaches the ground
    end
end

x_ft = x*m2ft;
y_ft = y*m2ft;

time_of_flight_s            %units in s
max_height_ft = max(y_ft)   %units in ft
range_ft = xf*m2ft          %units in ft
final_speed_mph = vf/mph2mps    %units in mph
energyLost_J = 0.5*m*(v0^2 - vf^2)  %units in J

%the time taken, range, and maximum height of the ball are the same as what
%was found in the excel file

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
title({'ECE 202, Project 2, Phase 5: Trajectory of a baseball', ...
    'drag vs without drag'}, 'FontSize', 22)
xlabel('x (ft)', 'FontSize', 18)   
ylabel('y (ft)', 'FontSize', 18)
str1 = sprintf('with drag and C = %g',C);
legend({'without drag', str1}, ...
    'FontSize', 18)

T0 = 5.7;
H0 = 114;
R0 = 446;

error_percentage_time_of_flight = 100*(time_of_flight_s - T0)/T0
error_percentage_max_height = 100*(max_height_ft - H0)/H0
error_percentage_range = 100*(range_ft - R0)/R0

%The error percentages for maximum height and time of flight are higher
%than the range of the baseball. Which means this value of C is only
%useful to estimate the range for the baseball. A scope of error can be
%that in real life the ball was struck from a height greater than 0, which
%was the starting point taken in this investigation.


