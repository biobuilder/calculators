% calculates ra given some periapsis circular orbit value and a multiplier
% used for positioning relay satellites
% Damon Printz
% 2/12/2021

function Ra = space_circular_orbit (Rp, k)
Ra = Rp * (2*k^(2/3)-1);
endfunction
