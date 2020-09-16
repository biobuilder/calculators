% Centrifugal pump fundamental calculations
% Used to calculate specific pump types 
% Damon Printz
% 8-25-2020

close all; clc; clear;

% input pump characteristics

rv = 5000 * 2*pi()/60; % operating rotational velocity (rad/s)
Q = 10E-6; % fluid flow rate (m^3/s)
h = 2; % head (m)
W = 25; % power (watts)
rho = 1.225 * 1000; % fluid density (kg/m^3)
g = 9.81; % gravity (m/s^2)

beta2 = pi()/2; % outlet blade angle (rad)


% calculate outer diameter of radial impeller pump
r2 = sqrt(h*g/rv^2);
fprintf('outer radius = %f m\n',r2)

% calculate tangential tip velocity 
U2 = r2 * rv % (m/s)

% calculate power required to achieve the given flow rate 
W = rho*Q*U2^2

% calculate width of impeller blades at outlet



% calculate inlet hole diameter 


















