function dist = distance(pts)
  # finds distance between pairs of points and adds them sequentially
  # always absolute value of distance
  # must be at least two points
  # pts rows correspond to the different axes
  dist = 0;
  
  # iterate from one point to the next
  for i = 1:size(pts,1)-1
    # distance is the square root of the coordinate differences squared and summed together
    dist = dist + sqrt(sum((pts(i+1,:)-pts(i,:)).^2));
  endfor
  
  
  
  
endfunction
