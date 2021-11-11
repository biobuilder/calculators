% time calculation script
% determines orbital radii for orbits of specified period difference
% used for positioning constellations of satellites


r = 4000000; % initial circular radius r 
ratio = 35.3553391/100; % ratio of T1 to T2

coefs = [ (1),                % third degree constant
          (2+r),              % second degree constant
          (2*r+r^2),          % 1st degree constant 
          (1-8*(ratio)^2)*r^3]; % zeroth degree constant 

roots(coefs)