## Copyright (C) 2020 damon
## 
## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see
## <https://www.gnu.org/licenses/>.

## -*- texinfo -*- 
## @deftypefn {} {@var{retval} =} orbital_plotter (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: damon <damon@FUNRIG>
## Created: 2020-04-02


% 

function [retval] = orbital_plotter (mu, Rp, Ra, inc, argp, lonasc)


% Orbital plotter
% Damon Printz
% 3/14/2020

% Plots orbits of planets in the solar intertial frame 
% translates a given radius and true anomaly into X,Y,Z coordinates in the intertial frame 


% input orbital parameters
mu = mu_Sun;
Rp = Rp_Ker; % periapsis 
Ra = Ra_Ker; % apoapsis
inc = inc_Ker; % inclination (rad)
argp = argp_Ker; % argument of periapsis (rad)
lonasc = lonasc_Ker; % longitude of ascending node (rad)

% get orbital parameters from this definition
[Vp,Va,T,e,p,E,h] = elliptical_param (mu, Rp, Ra)

% generate a list of true anomaly angles 
numPts = 200;
tanom = linspace(0, 2*pi(), numPts);
% generate radii lengths as a function of these angles 
rlen = p./(1+e.*cos(tanom));

% obtain orbital circle radii in inertial frame 
% convert to x, y, z coordinates in the inertial frame 
r = zeros(length(tanom),3);

r(:,1) = (cos(lonasc).*cos(argp+tanom(:))-sin(lonasc).*cos(inc).*sin(argp+tanom(:)));
r(:,2) = (cos(lonasc).*cos(inc).*sin(argp+tanom(:))+sin(lonasc).*cos(argp+tanom(:)));
r(:,3) = (sin(inc).*sin(argp+tanom(:)));

% calculate orbital radius vectors
r_orb = zeros(size(r));
for(i = 1:1:length(rlen))
r_orb(i,:) = r(i,:) .* rlen(i);
endfor 

% calculate the ascending node vector 
ra = zeros(3,1);
ra(1) = cos(lonasc);
ra(2) = sin(lonasc);
rd = -1.*ra; % descending node vector 

% find the index of the orbital point closes to this vector 
minRa = 1;
minRaIndex = 1;
minRd = 1;
minRdIndex = 1;
for(i = 1:1:length(tanom))
  distRa = sqrt((r(i,1)-ra(1))^2+(r(i,2)-ra(2))^2+(r(i,3)-ra(3))^2);
  distRd = sqrt((r(i,1)-rd(1))^2+(r(i,2)-rd(2))^2+(r(i,3)-rd(3))^2);
  if(distRa < minRa)
    minRaIndex = i;
    minRa = distRa;
  endif 
  if(distRd < minRd)
    minRdIndex = i;
    minRd = distRd;
  endif 
endfor 

% plot the orbital radius vectors
hold on 
plot3(r_orb(:,1),r_orb(:,2),r_orb(:,3))

% plot the apoapsis, periapsis, ascending node, etc.
hold on 
% periapsis 
plot3(r_orb(1,1),r_orb(1,2),r_orb(1,3), '*')
% apoapsis
plot3(r_orb(floor(numPts/2),1),r_orb(floor(numPts/2),2),r_orb(floor(numPts/2),3), '*')
% ascending node 
plot3(r_orb(minRaIndex,1),r_orb(minRaIndex,2),r_orb(minRaIndex,3),'*')
% descending node
plot3(r_orb(minRdIndex,1),r_orb(minRdIndex,2),r_orb(minRdIndex,3),'*')
legend('orbit','apoapsis','periapsis','ascending node','descending node')

retval = 1;

endfunction
