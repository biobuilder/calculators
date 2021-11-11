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
## @deftypefn {} {@var{retval} =} get_vel (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: damon <damon@FUNRIG>
## Created: 2020-03-06

% input mu, Ra, Rp and desired orbital intersect radius 
% get velocity vector as polar coordinates radial out, tangent
% get true anomaly


function [V, ta] = get_vel (mu, Rp, Ra, R)

  V = zeros(1,2); % radial out, perpendicular

  % get orbital elements from the elliptical orbit 
  [Vp,Va,T,e,p,E,h] = elliptical_param(mu, Rp, Ra);
  a = (Ra+Rp)/2;
  
  % true anomaly 
  ta = acos((a*(1-e^2)/R-1)/e);
  
  % get speed at the desired radius 
  speed = sqrt(2*(E+mu/R));
  
  % calculate angle between local horizontal and velocity vector 
  localAng = acos(h/(R*speed));
  
  % break into vector form in polar coordinates 
  V(1) = -speed*sin(localAng); % radial out 
  V(2) = speed*cos(localAng); % perpendicular 


endfunction
