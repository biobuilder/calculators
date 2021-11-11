% builds radii vectors from orbital elements
% Damon Printz
% 4/5/2020

% input all orbital elements
% must include true anomaly 


% inputs 
% tanom = true anomaly [rad] 
% inc = inclination [rad] 
% argp = argument of periapsis [rad]
% lonasc = longitude of ascending node [rad]
% e = eccentricity
% p = semilatus rectum

% Outputs 
% r     = orbital vectors 
% ra    = ascending node vector 
% rd    = descending node vector    
function [r, ra, rd] = build_rvect_from_elements (tanom, inc, argp, lonasc, e, p)

  % obtain orbital radii in inertial frame 
  % convert to x, y, z coordinates in the inertial frame 
  
  % calculate the length of radii at the given true anomaly angles
  rlen = p./(1+e.*cos(tanom));
  
  % initialize an orbital radii array 
  r = zeros(length(tanom),3);

  % convert true anomaly coordinates into the inertial reference frame
  % x, y and z coordinates
  r(:,1) = (cos(lonasc).*cos(argp+tanom(:))-sin(lonasc).*cos(inc).*sin(argp+tanom(:))) .* rlen(:);
  r(:,2) = (cos(lonasc).*cos(inc).*sin(argp+tanom(:))+sin(lonasc).*cos(argp+tanom(:))) .* rlen(:);
  r(:,3) = (sin(inc).*sin(argp+tanom(:))) .* rlen(:);

  
  % calculate the ascending node unit vector 
  ra = zeros(3,1);
  ra(1) = cos(lonasc);
  ra(2) = sin(lonasc);
  % calculate the descending node unit vector 
  rd = -1.*ra;
  
endfunction
