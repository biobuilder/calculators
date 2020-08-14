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
## @deftypefn {} {@var{retval} =} q_calc (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: damon <damon@FUNRIG>
## Created: 2020-03-06

% dynamic pressure assist calculator
% select a planet, input altitude Alt (m) and speed vel(m/s)
% outputs dynamic pressure q (Pa/m^2), mach number M and speed of sound a (m/s)

function [q,M,a] = q_calc(alt,vel,planet)
  
  % get atmospheric conditions on the planet
  [T,p,rho,gam,R] = atmosphere(alt,planet);
  % get speed of sound 
  a = sqrt(gam*R.*T);
  % get mach number 
  M = vel./a;
  % get dynamic pressure 
  q = 0.5*gam.*p.*M.^2;
  
endfunction
