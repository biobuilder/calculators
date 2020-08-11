function [center] = centroid(pts)
  % calculates the centroid of an arbitrary number of points
  % pts is a 1 to 3 - column matrix, where each row represents a new coordinate
  % column 1 = x, column 2 = y, column 3 = z
  center = zeros(1,size(pts,2)); 
  
  % basically, all you do is average the coordinate values
  for i = 1:size(pts,2)
      center(1,i) = sum(pts(:,i)) / size(pts,1);
  endfor
  
endfunction
