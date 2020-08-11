# Space Engineers coordinate generation scripts
# Damon Printz
# 2/6/2019

# These zones will help you calculate GPS coordinates for different solutions
clc; clear; close all;


# calculate midpoint between two coordinates
v1 = [16625.13, 145845.13, -114121.81]; 
v2 = [15837.70, 126965.45, -112604.22];

dist = norm(v1-v2);

fprintf('distance = %.2f km\n',dist/1000)

# point spacing
spacing = 18960/2; # distance between the points

# Output the points for spacing

numCoords = ceil(dist / spacing)+1;

% initialize an array containing all the coordinates
% starts at vector 1 and ends at vector 2
coords = zeros(numCoords,3);

coords(1,:) = v1; % first coordinate
coords(numCoords,:) = v2; % last coordinate


% calculate coordinate set inbetween
for i = [2:numCoords-1]
  coords(i,:) = coords(1,:) - ((i-1)*spacing).* ... 
                (coords(1,:)-coords(numCoords,:)) ./ ... 
                norm(coords(1,:) - coords(numCoords,:))
endfor


# print the new vector set
fprintf('x, y, z\n')
for i = [1:numCoords]
  fprintf('%.0f, %.0f, %.0f\n',coords(i,:))
endfor

