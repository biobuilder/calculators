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
## @deftypefn {} {@var{retval} =} altcalc (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: damon <damon@FUNRIG>
## Created: 2020-02-27

function [pressure, density] = altcalc(altitude, type)

% altitude (m), pressure (Pa), pressure (atm), density (kg/m^3)
altData = zeros(17,4);
maxAlt = 0; 

% altitude data for Kerbin
if(type == 'k')
altData = ...
[0	    101325	  1.000	    1.225; ...
2500	  69015	    0.681	    0.898; ...
5000	  45625	    0.450	    0.642; ...
7500	  29126	    0.287	    0.446; ...
10000	  17934	    0.177	    0.288; ...
15000	  6726	    0.066	    0.108; ...
20000	  2549	    0.025	    0.040; ...
25000	  993.6	    0.010	    0.015; ...
30000	  404.1	    0.004	    0.006; ...
40000	  79.77	    0.001	    0.001; ...
50000	  15.56	    0.000	    0.000; ...
60000	  2.387	    0.000	    0.000; ...
70000	  0.000	    0.000	    0.000];

maxAlt = 70000;

% Eve planetary body 
elseif type == 'e'
altData = [...
0	      506625	  5.000
2500	  403913	  3.986
5000	  316277	  3.121
7500	  242677	  2.395
10000	  182072	  1.797
12500	  133423	  1.317
15000	  95689	    0.944
20000	  44335	    0.438
25000	  18073	    0.178
30000	  8233	    0.081
35000	  5196	    0.051
40000	  3500	    0.035
50000	  121.8	    0.001
60000	  23.00	    0.000
70000	  4.344	    0.000
80000	  0.8205	  0.000
90000	  0	        0.000];

maxAlt = 90000;
% Duna planetary body 
elseif type == 'd'

% Jool planetary body 
elseif type == 'j'

% Laythe planetary body 
elseif type == 'l'


endif

% calculates the air pressure and density 
% at different altitudes for the given dataset
if(altitude < maxAlt && altitude >= 0)
  pressure = interp1(altData(:,1),altData(:,2),altitude);
  density = interp1(altData(:,1),altData(:,4),altitude);
else 
  pressure = 0;
  density = 0;
endif 

endfunction
