function [pts] = orbitize(center, axis, initial, radius, npts)
  % orbitize creates a series of points equally spaced around the given axis
  % center is the center of the rotation
  % axis is the coordinate representing the positive axis direction
  % initial is the first coordinate direction (projected onto a line perpendicular to the axis)
  % radius is the projection radius
  % npts is the number of points 
  
  % calculate right ascension of the ascending node
  raan = cross((axis-center), (initial-center));
  % calculate initial coordinate direction
  iniVect = cross(raan, (initial-center));
  
  
  
endfunction
