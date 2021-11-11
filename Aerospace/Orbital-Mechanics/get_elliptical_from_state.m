% calculates the new elliptical orbital elements from a single point measurement 
% outputs
% Damon Printz
% 4/5/2020
% Rp = periapsis, Ra = apoapsis, Vp = periapsis velocity, Va = apoapsis velocity 
% T = orbital period, e = eccentricity, p = semilatus rectum, E = specific energy
% inputs 
% mu = gravitational constant, R = current radius, V = current velocity,
% h = current 

function [Rp,Ra,Vp,Va,T,e,p,E] = get_elliptical_from_state(mu, R, V, h)

E = V^2/2-mu/R;

aq = -E;
bq = -mu;
cq = h^2/2;

if(bq^2-4*aq*cq < 0)
  fprintf('get_apoapsis_from_periapsis: incompatible orbital data\n')
else 
  Ra = (-bq+sqrt(bq^2-4*aq*cq))/(2*aq);
  Rp = (-bq-sqrt(bq^2-4*aq*cq))/(2*aq);
endif 

Va = h/Ra;
Vp = h/Rp;

a = (Ra+Rp)/2; 
T = 2*pi()*sqrt(a^3/mu);
e = (Ra-Rp)/(Ra+Rp);
p = a*(1-e^2);

endfunction
