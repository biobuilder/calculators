% outputs the velocity at the periapsis of a hyperbolic trajectory
% based on V_inf (hyperbolic excess velocity)
% Damon Printz

function [Vp] = hyperbolic_periapsis (mu, Rp, Vinf)

Vp = sqrt(Vinf^2+2*mu/Rp);

endfunction
