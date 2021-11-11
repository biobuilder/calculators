% Rocket Launch System
% Purpose: develop a gravity turn schedule
% Damon Printz
% 2/27/2020
clc; clear; close all;

% simulation constants 
tstep = 0.1; % time step in seconds 
startingLatitude = 0; % launch latitude in degrees
startingAltitude = 100; % launch altitude from sea level
startingVelocity = [100, 0]; % launch starting velocity (radial, horizontal)
startingAzimuth = 90; % launch azimuth on compass (degrees)
startingTip = 0; % starting attitude orientation (degrees)

maxSimSteps = 10000; % maximum simulation steps 

hold_attitude_for = [80,120/tstep]; % holds the thrust angle (1) [deg] at a given input for (2) steps 

% environment constants
Rb = 600000; % mean radius of body (m)
mu = 3.5316E12; % gravitational constant (m^3/s^2)
T = 21549.3; % sidereal day (s)
g0 = 9.81; % specific impulse gravitational constant (m/s^2)
type = 'k'; % planet type for altcalc

% payload
mPayload = 840+100+300+40; % mass of the payload (kg)
aPayload = pi()*1.25^2; % frontal area of the payload
cdRocket = 0.5; % coefficient of drag for the rocket

% staging information

fullMass = [2250*2+500, 1*7650]; % full mass of the stage
emptyMass = [25*2+500, 1*1500]; % empty mass of the stage (kg)
burnRate = [26.625, 1*145.65]; % fuel burn rate (kg/s)
isp = [345, 175]; % specific impulse (s)
areaStage = [0, 0*pi()*1.25^2]; % additional area from stage (m^2)

%{
fullMass = [1125*4+1500, 2*(1125*4+1200)]; % full mass of the stage
emptyMass = [125*4+1500, 2*(125*4+1200)]; % empty mass of the stage (kg)
burnRate = [68, 2*(78.9)]; % fuel burn rate (kg/s)
isp = [300, 265]; % specific impulse (s)
areaStage = [0, 3*pi()*1.32^2]; % additional area from stage (m^2)

fullMass = [1125*3+1500, 1*(1125*3+1200), 2*3560]; % full mass of the stage
emptyMass = [125*3+1500, 1*(125*3+1200), 2*750]; % empty mass of the stage (kg)
burnRate = [68, 1*(78.9), 2*119]; % fuel burn rate (kg/s)
isp = [250, 265, 170]; % specific impulse (s)
areaStage = [0, 0, 2*pi()*1.32^2]; % additional area from stage (m^2)

fullMass = [1125*4+1500, 2*(1125*4+1200)]; %3560*2]; % full mass of the stage
emptyMass = [125*4+1500, 2*(125*4+1200)]; %750*2]; % empty mass of the stage (kg)
burnRate = [68, 2*(78.9)]; % 119*2]; % fuel burn rate (kg/s)
isp = [250, 265]; % 170]; % specific impulse (s)
areaStage = [0, 2*pi()*1.32^2]; % additional area from stage (m^2)
%}

% simulation variables and calculations
mRocket = mPayload .* ones(maxSimSteps,1);
area = aPayload .* ones(maxSimSteps,1);
altitude = zeros(maxSimSteps,1); % spacecraft altitude from the mean radius of the planet
attitude = zeros(maxSimSteps,1); % spacecraft attitude relative to the surface of the planet

position = zeros(maxSimSteps,2); % radius, angle from starting vector(r, phi)
velocity = zeros(maxSimSteps,2); % velocity in radial and tangential directions 
%altitude(1) = startingAltitude; % starting from ground level
attitude(1) = startingTip * pi() / 180; % initial attitude
position(1,:) = [Rb + startingAltitude, 0]; % starting position in radial coordinates
velocity(1,:) = startingVelocity;
% forces 
thrust = zeros(maxSimSteps,1); % always inline with the rocket
drag = zeros(maxSimSteps,1); % always inline with the velocity vector 
weight = zeros(maxSimSteps,1); % always inline with the radial direction

% stages and mass 
stageCnt = length(fullMass); % current stage count (starts at the end and counts down to zero)
numStages = length(fullMass); % number of stages
curMass = fullMass; % current mass of each stage (kg)

endStep = maxSimSteps;

% calculate 10% marker for steps 
marker = floor(maxSimSteps / 10);

% run the simulation
% rules: every time a tank runs out of fuel, drop the stage and ignite the next 

for step = [1:1:maxSimSteps-1]
  
  % calculate the mass and frontal surface area of the rocket
  %mRocket(step) = mPayload;
  %area(step) = aPayload;
  
  % count down from current stage count to zero
  % add the masses and the areas of the stages 
  if(stageCnt >= 1)
      for stage = [stageCnt:-1:1]
          % if the current stage mass is still greater than the empty mass,
          % add its current mass to the total mass of the rocket 
          % otherwise decrease the stage count and do not add the mass
          if(curMass(stage) > emptyMass(stage))
            mRocket(step) = mRocket(step) + curMass(stage);
            area(step) = area(step) + areaStage(stage);
          else
            stageCnt = stageCnt - 1;
          endif
      endfor 
  endif
  
  
  % calculate force vectors
  % if the stage count is greater than the last stage,
  % burn some fuel and accelerate the rocket 
  if(stageCnt > 0)
    curMass(stageCnt) = curMass(stageCnt) - burnRate(stageCnt) * tstep;
    thrust(step) = g0 * isp(stageCnt) * burnRate(stageCnt);
  endif 
  
  % drag (always inline with velocity vector)
  if(position(step,1) >= Rb)
    [pressure, density] = altcalc(position(step,1)-Rb,type);
    drag(step) = -0.5 * area(step) * cdRocket * density * (velocity(step,1)^2 + velocity(step,2)^2);
  else
    endStep = step;
    disp('Rocket collided with the planet')
    break 
  endif
  
  % gravitational forces, aligned with radial out
  weight(step) = -mu * mRocket(step) * position(step,1)^(-2);

  % sum forces in r and p directions
  % attitude is held by thrusters for the given amount of steps 
  % otherwise the rocket is allowed to naturally curve over on its own 
  if(hold_attitude_for(2) > step)
    rForce = thrust(step)*sind(hold_attitude_for(1)) + drag(step)*sin(attitude(step)) + weight(step);
    pForce = thrust(step)*cosd(hold_attitude_for(1)) + drag(step)*cos(attitude(step));
  else 
    rForce = (thrust(step) + drag(step)) * sin(attitude(step)) + weight(step);
    pForce = (thrust(step) + drag(step)) * cos(attitude(step));
  endif
  
  % calculate change in velocity 
  % F = ma, a*t = DV
  rDV = (rForce / mRocket(step)) * tstep;
  pDV = (pForce / mRocket(step)) * tstep;
  velocity(step,:) = velocity(step,:) + [rDV, pDV];
  
  % calculate change in position 
  % V*t = DX
  drr = velocity(step,1) * tstep; % change in radial direction 
  dpp = velocity(step,2) * tstep; % change in tangential direction 

  newR = norm([position(step,1) + drr, dpp]);
  dPhi = atan(dpp / (position(step,1) + drr));
  newPhi = position(step,2) + dPhi;
  
  % calculate new position in radial coordinates 
  position(step+1,:) = [newR, newPhi];
  
  % rotate velocity vector to match new coordinates  
  velocity(step+1,1) = velocity(step,1)*cos(dPhi) + velocity(step,2)*sin(dPhi);
  velocity(step+1,2) = -velocity(step,1)*sin(dPhi) + velocity(step,2)*cos(dPhi);
  
  % calculate new attitude 
  % based on velocity vector
  if(velocity(step,2) > 0) 
    attitude(step+1) = atan(velocity(step,1) / velocity(step,2));
  elseif(velocity(step,1) > 0)
    attitude(step+1) = pi()/2; % 90 degrees straight up
  elseif(velocity(step,1) < 0)
    attitude(step+1) = -pi()/2; % 90 degrees straight down 
  else  
    attitude(step+1) = attitude(step); % no velocity, preserve prior attitude 
  endif

  
  % display progress status 
  if(mod(step,marker) == 0 || step == 1)
      percentDone = floor(100*step/maxSimSteps);
      fprintf('%i%%\n',percentDone)
  endif 
endfor

% clean up the data (get rid of the crash data) 
for step = [1:1:maxSimSteps]
    if(position(step,1) < Rb)
        position(step,1) = Rb; 
    endif 
endfor 

% create a time variable
time = [0:tstep:tstep*(endStep-1)];

disp('100% - finished!')

% Plots

%plot(position(1:endStep,2).*180/pi(), position(1:endStep,1))
%title('radius (m) vs. phi (deg)')

figure 
plot(time, position(1:endStep,1))
title('radius (m) vs. time (s)')
grid on

figure 
%plot(position(1:endStep,2).*180/pi(), sqrt(velocity(1:endStep,1).^2+velocity(1:endStep,2).^2))
%title('Speed (m/s) vs. position (deg)')
plot(time, sqrt(velocity(1:endStep,1).^2+velocity(1:endStep,2).^2))
title('speed (m/s) vs. time (s)')
grid on

figure
plot(time, thrust(1:endStep))
title('thrust (N) vs. time (s)')
grid on

figure
plot(time, drag(1:endStep))
title('drag (N) vs. time (s)')
grid on

figure 
%plot(position(1:endStep,2).*180/pi(), mRocket(1:endStep))
%title('Rocket mass (kg) vs. position (deg)')
plot(time, mRocket(1:endStep))
title('rocket mass (kg) vs. time (s)')
grid on

figure
plot(time, attitude(1:endStep).*180/pi())
title('attitude (deg) vs. time (s)')
grid on