% calculates minimum (ideal) Delta V to get to orbit 
% from given gravitational constant mu, starting radius Re,
% and orbital radius Ro
% does not include velocities due to the rotation of the body
% Damon Printz

function DV = to_orbit (mu, Re, Ro)
Eo = -mu/(2*Ro);
Ee = -mu/Re;
DE = (Eo-Ee);
DV = sqrt(2*DE);
endfunction
