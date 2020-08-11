function [gps] = encode(coord, name)
  # encodes a point into a GPS coordinate string for Space Engineers
  # input coordinates coord, where columns represent x, y, z
  # input a name for the GPS coordinate
  # output a Space Engineers GPS String gps
  # Example string: 'GPS:E - Earth Lake Lavon:-18368.01:19464:54731.27:'
  
  gps = ['GPS:',name,':',num2str(coord(1,1)),':',num2str(coord(1,2)),':',num2str(coord(1,3)),':'];
  
endfunction
