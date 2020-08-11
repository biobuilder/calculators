function [coord] = decode(gps)
  # decodes a Space Engineers GPS string and converts it to a coordinate
  # input GPS string gps
  # get out a coordinate, where columns represent x,y,z
  # Example string: 'GPS:E - Earth Lake Lavon:-18368.01:19464:54731.27:'
  
  # initialize coordinate variable
  coord = zeros(1,3);
  # split string into constituent parts 
  pieces = strsplit(gps,':');
  # convert the required parts into coordinates
  coord(1,1) = str2num(char(pieces(1,3)));
  coord(1,2) = str2num(char(pieces(1,4)));
  coord(1,3) = str2num(char(pieces(1,5)));
  
endfunction
