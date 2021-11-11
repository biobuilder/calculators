%
%
%                          Launch Windows 
%
%       Generates a list of launch windows between two planets
%   Goals: 
%   1. To determine optimal launch and return windows 
%      for interplanetary trajectories
%      
%   2. To determine the trajectory characteristics
%
%   Created by:     Damon Printz
%   Started on:     4/5/2020
%   Current update: 4/5/2020
% 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
clear; clc; close all;

% Must run 'ksp_constants' first
source('ksp_constants.m')

% Simulation Constants
muSun = mu_Sun; % gravitational parameter [m^3/s^2]

% orbital mechanics structures
host = struct("mu",     mu_Ker,...    % gravitational parameter [m^3/s^2]
              "ra",     Ra_Ker,...    % apoapsis [m]
              "rp",     Rp_Ker,...    % periapsis [m]
              "inc",    inc_Ker,...   % inclination [rad]
              "argp",   argp_Ker,...  % argument of periapsis [rad]
              "lonasc", lonasc_Ker,...% longitude of ascending node [rad]
              "torb",   Torb_Ker,...  % orbital period [s]
              "ma",     MA_Ker,...    % starting mean anomaly at 0s UT [rad]
              "e",      0, ...        % eccentricity (calculated later)
              "p",      0 ...         % semilatus rectum (calculated later) [m]
              );
              
targ = struct("mu",     mu_Moh,...    % gravitational parameter [m^3/s^2]
              "ra",     Ra_Moh,...    % apoapsis [m]
              "rp",     Rp_Moh,...    % periapsis [m]
              "inc",    inc_Moh,...   % inclination [rad]
              "argp",   argp_Moh,...  % argument of periapsis [rad]
              "lonasc", lonasc_Moh,...% longitude of ascending node [rad]
              "torb",   Torb_Moh,...  % orbital period [s]
              "ma",     MA_Moh,...    % starting mean anomaly at 0s UT [rad]
              "e",      0, ...        % eccentricity (calculated later)
              "p",      0 ...         % semilatus rectum (calculated later) [m]
              );

% Simulation constants - time
tstart = Torb_Ker * 1 + (284+86)*6*3600;   % starting UT time [s]
tend   = tstart + Torb_Ker * 2;         % ending UT time [s]
tstep  = 2*3600;                        % time step [s]

% Total trajectory-solver steps to perform 
% Simulated positions for planets must be calculated ahead of time. 
% The target will not have valid trajectories near the end of its time 
% iterable.  Therefore the total number of steps may be less than 
% the total number of time steps.  
totalSteps = floor( ( ( tend - tstart ) / tstep ) / 2 )

% Opposition angles 
% Determine whether or not a trajectory will be 
% attempted.
% The angle between the departure and arrival vectors 
% must not exceed either boundary.
% (The angle will be configured between 0 and 2*pi()) 
minOppAng = (180 - 15) * pi()/180;  % [rad]
maxOppAng = (180 + 15) * pi()/180;  % [rad]

% Arrival time multipliers
% Determine whether or not a trajectory will be
% attempted.
% The target must arrive between the boundary values
% determined by multiplying the min and max constants 
% by 1/2 of a standard Hohmann transfer time period 
% with semimajor axis equal to the average length of 
% the departure and arrival vectors.
minArriveTimeMult = 0.5;
maxArriveTimeMult = 1.5;

% Inclination boundary 
% Determine whether or not a trajectory will be 
% attempted.
% The inclination of the cross product of the departure  
% and arrival vectors must be less than this amount 
% in order to be considered valid.
% This is because high inclination orbits are difficult
% to achieve and are not worth including in the analysis.
minInc = -30 * pi()/180;    % [rad]
maxInc =  30 * pi()/180;    % [rad]

% Trajectory analysis constants 
% The trajectory analysis must obtain a trajectory with 
% a transit time within timeAccu of the target transit 
% time to be considered valid  
timeAccu = 3600; % [s]
% maximum number of attempts to achieve a valid trajectory
maxTimeIter = 15;
% maximum Delta V allowed for a trajectory 
maxDV = 100000; % [m/s]

% True anomaly newton solver constants
% minimum step accuracy 
% The solver exits when step size is less than this value.
tanomAccu = 0.01 * pi()/180; % [rad]
% maximum steps before exiting 
tanomMaxSteps = 15;


% create variables storing trajectory information
% current trajectory index count 
trajIndCnt      = 0;
% interplanetary transfer window time indexes
hostTimeInd     = zeros( 1, totalSteps ); % departure time vectors 
trajTimeInd     = zeros( 1, totalSteps ); % arrival time vectors
% hyperbolic excess velocity at departure 
trajHypDep      = zeros( 1, totalSteps );
% hyperbolic excess velocity at arrival 
trajHypArr      = zeros( 1, totalSteps );
% inclination of interplanetary trajectories 
trajinc         = zeros( 1, totalSteps );
% longitude of ascending node of interplanetary trajectories 
trajlonasc      = zeros( 1, totalSteps );
% eccentricity of interplanetary trajectory
traje           = zeros( 1, totalSteps );
% launch point of interplanetary trajectories 
trajlp          = zeros( 1, totalSteps );


% - - - - - - - - - - - - - Start simulation - - - - - - - - - - - - - - - - - -

% fill out the rest of the orbital constants for the bodies
[Vp,Va,T,host.e,host.p,E,h] = elliptical_param (muSun, ...
                                                host.rp, ...
                                                host.ra );
                                                
[Vp,Va,T,targ.e,targ.p,E,h] = elliptical_param (muSun, ...
                                                targ.rp, ...
                                                targ.ra );

% 1. Generate a set of points in space representing the host
%    and target body vectors from the inertial reference frame
%    centered on the sun.

time = [tstart:tstep:tend]; % time values

% generate true anomaly angles from these time steps for each
% body 
tanomHost = zeros(size(time));
tanomTarg = zeros(size(time));

% use a newton solver to get true anomaly vs. time 
for( i = [1:1:length(time)])

  % - - - - - - - Host - - - - - - - -
  % calculate mean anomaly at time(i)
  % mean anomaly progression at time t + initial mean anomaly
  meanAnomHost = (2*pi()/host.torb) * time(i) + host.ma;
  % convert to a number between 0 and 2*pi()
  while(meanAnomHost > 2*pi())
    meanAnomHost = meanAnomHost - 2*pi();
  endwhile 
  % use newton solver to get eccentric anomaly at time(i)
  eccAnomHost = newton_Ecc_from_Mean( meanAnomHost, ... % mean anomaly [rad]
                                      host.e, ...       % eccentricity
                                      tanomAccu, ...    % accuracy [rad]
                                      tanomMaxSteps ... % max steps
                                      );
                                      
  % calculate and store true anomaly in the array 
  tanomHost(i) = acos( ( cos( eccAnomHost ) - host.e ) / ...
                       ( 1 - host.e * cos( eccAnomHost ) ) );
  
  % run backwards if mean anomaly >= pi()
  if(meanAnomHost >= pi())
    tanomHost(i) = 2*pi() - tanomHost(i);
  endif 
  
  % - - - - - - - Target - - - - - - -
  % calculate mean anomaly at time(i)
  % mean anomaly progression at time t + initial mean anomaly
  meanAnomTarg = (2*pi()/targ.torb) * time(i) + targ.ma;
  % convert back to a number between 0 and 2*pi()
  while(meanAnomTarg > 2*pi())
    meanAnomTarg = meanAnomTarg - 2*pi();
  endwhile 
  % use newton solver to get eccentric anomaly at time(i)
  eccAnomTarg = newton_Ecc_from_Mean( meanAnomTarg, ... % mean anomaly [rad]
                                      targ.e, ...       % eccentricity
                                      tanomAccu, ...    % accuracy [rad]
                                      tanomMaxSteps ... % max steps
                                      );
                                      
  % calculate and store true anomaly in the array 
  tanomTarg(i) = acos( ( cos( eccAnomTarg ) - targ.e) / ...
                       ( 1 - targ.e * cos( eccAnomTarg ) ) );
  
  % run backwards if mean anomaly >= pi()
  if(meanAnomTarg >= pi())
    tanomTarg(i) = 2*pi() - tanomTarg(i);
  endif   
  
endfor 


% generate orbital vectors
% r = orbital vectors
[ rHost ] = ...
  build_rvect_from_elements (tanomHost, ...  % true anomalies [rad]
                             host.inc, ...   % inclination [rad]
                             host.argp, ...  % argument of periapsis [rad]
                             host.lonasc, ...% longitude of ascending node [rad]
                             host.e, ...     % eccentricity
                             host.p ...      % semilatus rectum [m]
                             );
                                       
[ rTarg ] =  ...
  build_rvect_from_elements (tanomTarg, ...  % true anomalies [rad]
                             targ.inc, ...   % inclination [rad]
                             targ.argp, ...  % argument of periapsis [rad]
                             targ.lonasc, ...% longitude of ascending node [rad]
                             targ.e, ...     % eccentricity
                             targ.p ...      % semilatus rectum [m]
                             );


                                 
%  - - - - - - - - - - - - Run Trajectory Analysis - - - - - - - - - - - - - - -                             

fprintf('radii for planets calculated\n')

% start at the 2nd iterable for velocity estimation (centered on the point)
for(i = [2:1:totalSteps])

  % store the current position vector 
  r1 = rHost(i,:);
  r1len = norm(r1);
  
  % generate a list of potential target vectors 
  
  fprintf('searching for target vectors, iteration %i\n',i)
  
  % create a list of indexes from which to store target indexes of interest  
  vind = zeros(1,totalSteps); % array of valid target indexes
  vindi = 0; % current valid index index count (0 means did not pass)
  
  % generate a list of target vectors 
  % first find all targets between boundary angles
  for(g = [2:1:length(time)])
      
      % store the potential target vector 
      r2 = rTarg(g,:);
      r2len = norm(r2);
      
      % calculate the angle between the two vectors 
      angle = acos(dot(r1, r2) / (r1len * r2len));
      
      if(angle < 0)
        angle = angle + 2*pi();
      endif
      
      % check to see if the angle is valid
      % if so, proceed with the next test
      % otherwise, proceed to the next test point - inclination check
      if(angle > minOppAng && angle < maxOppAng)
        
        % check to see if the inclination is valid 
        % calculate the normal vector of the orbit 
        testNormal = cross(r1,r2);
        % determine the inclination 
        testInc = acos( dot(testNormal, [0,0,1]) / norm(testNormal) );
        
        % if the inclination is valid, proceed with the next test 
        % transit time must be within a certain range to be considered viable 
        if(testInc > minInc && testInc < maxInc)
          
          % make sure transit time is within bounds 
          testHohmannTime = pi()*sqrt(((r1len+r2len)/2)^3/muSun);
          % generate time boundaries for the arrival of the object 
          % must be the departure time at iteration i
          % plus the hohmann transfer time multiplied by the boundaries 
          minTime = time(i) + testHohmannTime * minArriveTimeMult;
          maxTime = time(i) + testHohmannTime * maxArriveTimeMult;
          
          % the time at index (g) must be within bounds 
          if(time(g) > minTime && time(g) < maxTime)
            
            % if so, all checks have passed
            % this index will be stored as a potential workload
            % for the current assessment
            vindi = vindi + 1;
            vind(vindi) = g;
            
          endif
          
        endif 
        
      endif 
  
  endfor

  % begin trajectory analysis for this job if valid target indeces were
  % discovered
  if(vindi > 0)
    
    % output number of targets found 
    fprintf('%i targets\n', vindi)
    
    % create variables to store the generated trajectory data 
    % for later reference 
    % hyperbolic excess velocity at departure [m/s]
    hexcDep = zeros(1, vindi);
    % hyperbolic excess velocity at arrival [m/s]
    hexcArr = zeros(1, vindi);
    % inclination of target trajectory [rad]
    inct = zeros(1, vindi);
    % argument of periapsis of all trajectories [rad]   
    argpt = zeros(1, vindi);
    % longitude of ascending node of all trajectories [rad] 
    lonasct = zeros(1, vindi);
    % eccentricity of all trajectories 
    et = zeros(1, vindi);
    % launch point (ascending node)
    lpt = zeros(1, vindi);
    
    % iterate through the list of target indices of interest 
    for(f = [1:1:vindi])
      
      % grab the target vector using the index specified in the index array 
      % of interest 
      r2 = rTarg(vind(f),:);
      r2len = norm(r2);
      
      % calculate the velocity vector of the departure and arrival planets 
      % at this point through 2-point centered estimation 
      % average displacement / total time
      vp1 = (rHost(i+1, :) - rHost(i-1, :)) ./ (2*tstep);
      vp2 = (rTarg(vind(f)+1, :) - rTarg(vind(f)-1, :)) ./ (2*tstep);
      
      % determine the least expensive trajectory for the given time
      %  1. determine DV escape and DV capture 
      %  2. find minimum total DV
      %  3. log the best valid trajectory information
      
      % transit time (time until the target reaches the specified vector)
      % this is a measurement of the difference in time between departure 
      % and the time at which the target will reach the destination point 
      tt = time(vind(f)) - time(i);
     
      % get orbital elements of the trajectory 
      [hexcDep(f), ...    % departure hyperbolic excess velocity 
       hexcArr(f), ...    % arrival hyperbolic excess velocity 
       inct(f), ...       % trajectory inclination [rad]
       argpt(f), ...      % trajectory argument of periapsis [rad]
       lonasct(f), ...    % trajectory longitude of ascending node [rad]
       et(f), ...         % trajectory eccentricity
       lpt(f) ...         % launch point (ascending node) [rad]
       ] = solve_gaussian_trajectory (...
       r1, ...            % departure vector 
       r2, ...            % arrival vector 
       muSun, ...         % gravitational parameter
       tt, ...            % transit time 
       vp1, ...           % departure planet velocity 
       vp2, ...           % arrival planet velocity 
       timeAccu, ...      % trajectory time accuracy 
       maxTimeIter, ...   % maximum iterations
       maxDV ...          % maximum delta v allowed
       );       
      
    endfor
    
    % analyze the list of valid trajectories and select the ideal
    % value 
    % find minimum total DV 
    
    % ensure the trajectory is valid before logging it
    % the index of the best trajectory is default zero 
    bestIndex = 0;
    % this is the best DV value found thus far
    % set to be the highest value to ensure the smallest can be found  
    minDV = max(hexcDep(:) + hexcDep(:));
    
    % there must be a valid trajectory for a minimum to be found 
    if(minDV > 0)
      
      % iterate through the trajectory list and find the minimum 
      % hyperbolic excess velocity iteration
      % only accept nonzero excess velocities 
      for(f = [1:1:vindi])
        dvSum = hexcDep(f) + hexcArr(f);
        if(dvSum > 0 && dvSum < minDV)
          minDV = dvSum;
          bestIndex = f;
        endif
      endfor 
    
      % store the new trajectory information
      % if the best index was not incremented, then no good trajectories 
      % were calculated and this step must be skipped 
      if(bestIndex > 0)
        fprintf(' - valid trajectory found\n')
        % increment the current valid trajectory count 
        trajIndCnt = trajIndCnt + 1;
        % store the time index in the time index array
        % host time indexes  
        hostTimeInd(trajIndCnt) = i;
        % trajectory time indexes 
        trajTimeInd(trajIndCnt) = vind(bestIndex);
        % store the hyperbolic excess velocity at departure 
        trajHypDep(trajIndCnt) = hexcDep(bestIndex);
        % store the hyperbolic excess velocity at arrival 
        trajHypArr(trajIndCnt) = hexcArr(bestIndex);
        % store trajectory inclination 
        trajinc(trajIndCnt) = inct(bestIndex);
        % store trajectory argument of periapsis 
        trajargp(trajIndCnt) = argpt(bestIndex);
        % store trajectory longitude of ascending node 
        trajlonasc(trajIndCnt) = lonasct(bestIndex);
        % store trajectory eccentricity 
        traje(trajIndCnt) = et(bestIndex);
        % store trajectory launch point 
        trajlp(trajIndCnt) = lpt(bestIndex);
      endif
      
    endif 
    
  endif

endfor 


% output all the information 

% find best time of departure in the analysis period 
bstInd = 0;
bestVelocity = max(trajHypArr(:)+trajHypDep(:));
for( i = [1:1:trajIndCnt] )
  velSum = trajHypArr(i) + trajHypDep(i);
  if(velSum < bestVelocity)
    bestVelocity = velSum;
    bstInd = i;
  endif 
endfor 

% calculate best UT year, day, hour and minute 
bstTime = time(hostTimeInd(bstInd));
years = floor(bstTime / Torb_Ker);
days = floor((bstTime - years*Torb_Ker) / (6*3600));
hours = floor((bstTime - years*Torb_Ker - days*6*3600) / 3600);
minutes = floor((bstTime - years*Torb_Ker - days*6*3600 - hours*3600) / 60);
seconds = bstTime - years*Torb_Ker - days*6*3600 - hours*3600 - minutes*60;

newTime = [years, days, hours, minutes, seconds];
% optionally run Kerbin departure time for exact ascending node point 
if(host.torb == Torb_Ker )
  curTime = [years, days, hours, minutes, seconds];
  RA = trajlonasc(bstInd);
  [newTime, t_until, RAcur] = get_KSC_t_until_ra (RA, curTime);
endif 

if(bstInd > 0)
  % display best time of departure and transit time available for the analysis period 
  fprintf('\nBest index of departure = %i\n', bstInd )
  fprintf('UT time of departure = y%i, d%i, %ih:%im:%is\n', ...
           newTime(1), newTime(2), newTime(3), newTime(4), newTime(5) )
  fprintf('transit time = %.3f days\n', ...
           (time(trajTimeInd(bstInd)) - time(hostTimeInd(bstInd)))/(6*3600) )
  fprintf('Hyperbolic excess at departure = %.1f m/s, ', trajHypDep(bstInd) )
  fprintf('arrival = %.1f m/s, ', trajHypArr(bstInd) )
  fprintf('total = %.1f m/s\n', trajHypDep(bstInd) + trajHypArr(bstInd) )
  fprintf('Orbital elements:\n')
  fprintf('inclination = %.3f deg, ', trajinc(bstInd)*180/pi() )
  fprintf('longitude of ascending node = %.3f deg, ', trajlonasc(bstInd)*180/pi() )
  fprintf('argument of periapsis = %.3f deg, ', trajargp(bstInd)*180/pi() )
  fprintf('eccentricity = %.3f, ', traje(bstInd) )
  fprintf('launch point = %.3f deg\n', trajlp(bstInd)*180/pi() )
else 
  fprintf('No trajectories found\n')
endif 



% plot coordinate pairs
%{
figure
hold on 
plot3(0,0,0,'*','color','k')
for i = [1:1:trajIndCnt]
  plot3(rHost(hostTimeInd(i),1),rHost(hostTimeInd(i),2),rHost(hostTimeInd(i),3),'*','color','b')
  plot3(rTarg(trajTimeInd(i),1),rTarg(trajTimeInd(i),2),rTarg(trajTimeInd(i),3),'*','color','r')
endfor 
%}

% plot DV requirements vs. time 
figure 
plot( (time(hostTimeInd(1:trajIndCnt)) - tstart) ./ (6*3600), ...
      trajHypDep(1:trajIndCnt)+trajHypArr(1:trajIndCnt), '*')
grid on 
title('Total hyperbolic excess velocities vs. time')
xlabel('epoch time (days)')
ylabel('velocity [m/s]')

% plot DV requirements vs. time 
figure 
plot( (time(hostTimeInd(1:trajIndCnt)) - tstart) ./ (6*3600), ...
      trajinc(1:trajIndCnt)*180/pi(), '*')
grid on 
title('Departure inclination vs. time')
xlabel('epoch time (days)')
ylabel('inclination [deg]')

% plot transit times 
figure 
plot( (time(hostTimeInd(1:trajIndCnt)) - tstart) ./ (6*3600), ...
      (time(trajTimeInd(1:trajIndCnt)) - time(hostTimeInd(1:trajIndCnt)))/(6*3600), '*' )
title('Transit time vs. departure time')
xlabel('epoch time (days)')
ylabel('transit time (days)')
grid on 


% diagnostics 
%{
vtestp1 = zeros(totalSteps, 3);
vtestp2 = zeros(totalSteps, 3);

for(i = [2:1:totalSteps])
    vtestp1(i,:) = (rHost(i+1, :) - rHost(i-1, :)) ./ (2*tstep);
    vtestp2(i,:) = (rTarg(trajTimeInd(f)+1, :) - rTarg(trajTimeInd(f)-1, :)) ./ (2*tstep);
endfor 

figure 
plot(sqrt(vtestp1(2:totalSteps,1).^2+vtestp1(2:totalSteps,2).^2+vtestp1(2:totalSteps,3).^2))
hold on 
plot(sqrt(vtestp2(2:totalSteps,1).^2+vtestp2(2:totalSteps,2).^2+vtestp2(2:totalSteps,3).^2))

figure 
plot(sqrt(rHost(2:totalSteps,1).^2+rHost(2:totalSteps,2).^2+rHost(2:totalSteps,3).^2))

figure
plot3(rHost(1:totalSteps,1),rHost(1:totalSteps,2),rHost(1:totalSteps,3))
axis('square')
%}



























