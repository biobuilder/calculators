% Orbital Mechanics
% Basic Explorations
% Damon Printz
% 2/18/2020
clear; close all;
% Orbit Plotter for small satellites around a large body

% orbital constants
G = 6.67384E-11; % gravitational constant [m^3/(kg*s^2)]
ME = 5.972E24; % mass of the body [kg]
RE = 6378.1E3; % radius of the body [km]
mu = G*ME % gravitational parameter [m^3/s^2]

% orbital parameters
rp = 6500E3; % periapsis [m]
ra = 8000E3; % apoapsis [m]


% calculate orbital parameters
a = (rp+ra)/2 % semimajor axis [m]

e = (ra-rp)/(ra+rp) % eccentricity

p = rp*(1+e) % semi-latus rectum 

T = 2*pi()*sqrt(a^3/mu) % orbital period (s)

% calculate radii as a function of true anomaly
theta = linspace(0,2*pi(),100); % true anomaly [rad]
r = (a*(1-e^2))./(1+e.*cos(theta)); % radius [m]

% plot the orbit and points
polar(theta,r); % plot orbital path
hold on;
polar(0,rp,'x') % plot periapsis
polar(pi(),ra,'x') % plot apoapsis
polar(theta,ones(size(theta)).*RE) % plot host body 
title('orbital radius vs. true anomaly')
legend('orbit','periapsis','apoapsis','body')
