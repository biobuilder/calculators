% calculates escape velocity from a given radius and 
% standard gravitational parameter
% Damon Printz

function Vesc = escape (mu, R)
Vesc = sqrt(2*mu/R);
endfunction
