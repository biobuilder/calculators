% Centrifugal pump fundamental calculations
% Used to calculate specific pump types 
% Damon Printz
% 8-25-2020

close all; clc; clear;

% input pump characteristics

rv = 10000 * 2*pi()/60; % operating rotational velocity (rad/s)
Q = 10E-6; % fluid flow rate (m^3/s)
h = 10; % head (m)
%W = 25; % power (watts)
rho = 1.225 * 1000; % fluid density (kg/m^3)
g = 9.81; % gravity (m/s^2)

%beta2 = 40*pi()/180; % outlet blade angle (rad)


% calculate outer diameter of radial impeller pump
r2 = sqrt(h*(g-1)/(rv^2*(1-1/g)));
fprintf('outer radius = %f m\n',r2)

r2alt = sqrt(h*g/rv^2);
% calculate tangential tip velocity 
U2 = r2 * rv; % (m/s)
fprintf('blade tip velocity = %f m/s\n',U2)

% calculate power required to achieve the given flow rate 
Vo2 = U2; % assume tangential component is the same as blade velocity
W = rho*Q*U2*Vo2;
fprintf('power required = %f W\n',W)

% calculate width of impeller blades at outlet
Vr2 = U2*0.01; % assume exit velocity is small compared to blade velocity 
b2 = Q/(2*pi()*r2*Vr2); % (m)

fprintf('outlet blade width = %f\n',b2)
% calculate inlet hole diameter
Vr1 = Vr2*2; % ratio of inlet to outlet radial velocity component
b1 = b2;
r1 = Q/(2*pi()*b1*Vr1);
fprintf('inlet radius = %f m\n',r1)


















