% elliptical orbital parameters calculator
% input gravitational parameter, periapsis and apoapsis
% returns V periapsis, apoapsis, period, eccentricity, semilatus rectum,
% specific energy, specific momentum
% Damon Printz

function [Vp,Va,T,e,p,E,h] = elliptical_param (mu, Rp, Ra)
if(Ra == Rp)
  Va = sqrt(mu/Ra);
  Vp = Va;
else 
  Va = sqrt(2*mu*(1/Ra-1/Rp)/(1-Ra^2/Rp^2));
  Vp = Va*Ra/Rp;
endif 

a = (Rp+Ra)/2;
T = 2*pi()*sqrt(a^3/mu);
e = (Ra-Rp)/(Ra+Rp);
p = a*(1-e^2);
h = Rp*Vp;
E = Vp^2/2-mu/Rp;

endfunction
