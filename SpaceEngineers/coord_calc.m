# Space Engineers Main Coordinate Calculation Page
# Damon Printz
# 3/14/2019
clc; clear; close all;


%{
# Task: Find Lavon Entry Point, using a projection out to space
%gps1 = 'GPS:E - Earth Lake Lavon:-18368.01:19464:54731.27:' 
gps1 = 'GPS:C - The Dominion:-11472:12438:-58617:';
%gps_z = 'GPS:Z - Earth Lavon Zone:-30988.14:33713.83:94943.94:' # old entry zone coordinates
gps_z = 'GPS:C - Recycling Center:-840.59:28701.41:-97832.09:';
dropPoint = decode(gps1) # get actual drop point coordinates
ptz = decode(gps_z) # get old zone coordinates 
earthCenter = [0,0,0]; # center of the Earth coordinates

# calculate new entry zone point
#gravRadius = distance([earthCenter;ptz])
gravRadius = 107000; # gravitational radius of Earth
entryZone = project([earthCenter;dropPoint], gravRadius);

%newName = 'Z-Earth Lavon Entry Pt';
newName = 'Z - Dominion Entry Zone Zero G';
newGPSCoord = encode(entryZone, newName)
%}


%
gps2 = 'GPS:P - Mars Center:1029970:130274:1629717:'; # mars planet center
gps3 = 'GPS:P - Alien Planet Center:133707:126409:5733288:'; # alien planet center

pt1 = [0,0,0];
pt2 = decode(gps2);
pt3 = decode(gps3);

center = centroid([pt1;pt2;pt3])

gpsCenter = encode(center, 'C - Planetary Centroid')
%

%{
earthCenter = [0,0,0];
spawnPoint1 = 'GPS:S - Spawn Point 1:-12054.71:7397.01:-77821.58:';
gps_spawn = decode(spawnPoint1);
spawnRadius = distance([earthCenter; gps_spawn]); % get normal spawning radius
patrolRadius = 130000; % patrol satellite broadcast radius
# project out these coordinates into space

projects = [ 0,  0, -1;
            -1,  0, -1;
            -1,  0,  0;
            -1,  0,  1;
             0,  0,  1;
             1,  0,  1;
             1,  0,  0;
             1,  0, -1];

            
pts = zeros(size(projects));     
       
for i = 1:length(projects)
  pts(i,:) = project([earthCenter; projects(i,:)], patrolRadius);
  newName = ['C - Broadcast Zone ',num2str(i),'A'];
  newGPSCoord = encode(pts(i,:), newName);
  disp(newGPSCoord);
endfor

distance( [pts(1,:) ; pts(2,:)] )
%}
