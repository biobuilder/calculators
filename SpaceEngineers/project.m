function [coord] = project(pts, dist)
  # project from one coordinate through another coordinate a distance d
  # output the resulting coordinate
  # pts columns represent axes values
  # 1st row represents the starting coordinate
  # 2nd row represents the direction coordinate
  # dist represents the desired distance of the projection 
  
  coord = dist.*(pts(2,:)-pts(1,:))/sqrt(sum((pts(2,:)-pts(1,:)).^2)) + pts(1,:);
  
endfunction
