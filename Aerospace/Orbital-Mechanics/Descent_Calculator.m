% Descent calculator
% No Atmosphere
% Damon Printz
% 3/13/2020

clc;

mu = mu_Dun; %6.5138398E10; % gravitational parameter (m^3/s^2)
R_orbit = 60000+Req_Dun; % initial circular orbital radius (m)

r_mean = Req_Dun; % mean body radius (m)
body_rotation_period = Tsr_Dun; % rotational period of body (s)
landing_alt = 0; % landing altitude (m)
landing_lat = 0; % landing latitude (degrees)

deorbit_DE_ratio = .5; % initial deorbit burn to drop change in specific energy by the given ratio
touchdown_time = 30; % time in seconds until touchdown at landing_throttle setting 
landing_throttle = 1; % landing throttle setting (0 to 1 of peak throttle)
braking_throttle = 1; % braking throttle setting (0 to 1 of peak throttle)
deorbit_throttle = 1; % deorbiting throttle setting (0 to 1 of peak throttle)
Tacc = 1*9.81; % rocket acceleration at full thrust (m/s^2)

%g_kerb = 9.81; % local gravitational acceleration on Kerbin (m/s^2)
g_surface = mu/(r_mean+landing_alt)^2; % surface gravitational field constant (m/s^2)


% final descent stage (transition from point 3 to 4)
acc = (Tacc * landing_throttle - g_surface);
% make sure the acceleration actually slows down the rocket
if(acc > 0)
%V3f = (1-TWR*landing_throttle)*g_surface*touchdown_time; % descent velocity at point 3
%C = (TWR*landing_throttle-1)*g_surface*touchdown_time*(1-touchdown_time/2); % descent altitude (point 3 altitude)
V3f = acc*touchdown_time;
C = -(0.5*acc*touchdown_time^2-V3f*touchdown_time);

R3 = r_mean + landing_alt + C; % get radius of point 3

% initial deorbiting maneuver (transition from orbit to point 2)
E2i = -mu/(2*R_orbit); % specific energy of circular orbit 
E3f = V3f^2/2-mu/R3; % specific landing phase energy at point 3
DE2 = deorbit_DE_ratio * (E3f - E2i); % deorbiting change in specific energy 
V2i = sqrt(mu/R_orbit); % initial orbital velocity 
V2f = sqrt(2*(DE2+mu/(2*R_orbit))); % target velocity after reaching point 2
DV2 = V2f - V2i; % change in velocity required to deorbit the spacecraft


% braking maneuver at point 3 transitioning to slowly dropping to the surface 
% get energy of the new orbit 
E3i = V2f^2/2-mu/R_orbit; % new orbital energy in the 2->3 transition 
V3i = sqrt(2*(E3i+mu/R3)); % incoming velocity at point 3 
DV3 = V3f - V3i; % required Delta V change to enter landing phase 

% calculate change in true anomaly
[Rp,Ra,Vp,Va,T,e,p,E] = get_elliptical_from_state(mu, R_orbit, V2f, R_orbit*V2f);
a = (Ra+Rp)/2;
true_anomaly = acos((a*(1-e^2)/R3-1)/e);
delta_true_anomaly = (pi() - true_anomaly);

% calculate traversal time to surface 
ecc_anom = 2*atan(sqrt((1-e)/(1+e))*tan(true_anomaly/2)); % get eccentric anomaly of true anomaly
Mtraverse = ecc_anom - e*sin(ecc_anom); % get mean traversal anomaly
Ttraverse = T*Mtraverse/(2*pi()); % traversal time to get from true anomaly of 0 to current true anomaly from periapsis 
T23 = T/2 - Ttraverse; % get transition time from point 2 to 3

% estimate burn schedule/altitudes 
phi = acos(Ra*V2f/(R3*V3i)); % get local difference between tangent and orbital direction 
t_deorbit = abs(DV2)/(Tacc*deorbit_throttle);
t_brake = abs(DV3)/(Tacc*braking_throttle-g_surface);
braking_altitude = R3-r_mean-0.5*(Tacc*braking_throttle-g_surface)*sin(phi)*t_brake^2+V3i*sin(phi)*t_brake;

% output mission data 
fprintf('Burns\n')
fprintf('Deorbit DV = %.1f m/s,\nBraking DV = %.1f m/s,\nLanding DV = %.1f m/s,\nTotal DV = %.1f m/s\n',DV2,DV3,-V3f,DV2+DV3-V3f)

fprintf('\nVelocities\n')
fprintf('V in Orbit = %.1f m/s,\nV after deorbiting burn = %.1f m/s,\nV before brake = %.1f m/s,\nV after brake = %.1f m/s\n',V2i,V2f,V3i,V3f)

fprintf('\nAltitudes\n')
fprintf('Orbit = %.1f m,\nBraking burn = %.1f m,\nLanding Burn = %.1f m,\nSurface = %.1f m\n',R_orbit-r_mean,braking_altitude,R3-r_mean,landing_alt)

fprintf('\nBurn times\n')
fprintf('Deorbit burn time = %.1f s,\nBraking burn time = %.1f s,\nLanding burn time = %.1f s\n',t_deorbit,t_brake,touchdown_time)

fprintf('\nTime\n')
fprintf('End deorbit burn = %.1f s,\nStart brake burn = %.1f s,\nStart landing burn = %.1f s,\nTouchdown = %.1f s\n',t_deorbit,T23-t_brake,T23,T23+touchdown_time)

fprintf('\nTrue anomaly traversal over ground from initial burn = %.1f deg\n', delta_true_anomaly*180/pi())

else 
fprintf('thrust to weight ratio or throttle selection too low\n')
endif


