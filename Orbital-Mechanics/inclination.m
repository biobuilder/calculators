% calculates the Delta V required to change inclination of an orbit 
% Damon Printz
% input current velocity V and change in inclination di 
% outputs required burn vector [Prograde; Radial in; Normal]
% type controls what type of burn to make
% type = 0 is simply burn normal
% type = 1 is preserve all orbital parameters except inclination 

% reference: https://en.wikipedia.org/wiki/Orbital_inclination_change



function DV = inclination (V, di, type)

DV = zeros(3,1); % burn vector [Prograde; Radial in; Normal]

if(type == 0)
  DV(3) = V*tan(di); % only burn in normal direction
elseif(type == 1)
  DV(1) = V*(1-cos(di)); % retrograde burn (prograde as + direction)
  DV(3) = V*sin(di); % normal burn 
else 
  fprintf('inclination: invalid type\n')
endif 


endfunction
