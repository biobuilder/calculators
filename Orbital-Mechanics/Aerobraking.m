% Aerobraking calculator
% attempts to calculate DV when aerobraking on different planets 
% Damon Printz
% 3/6/2020

% outputs: DV change, reduction in apoapsis, peak G's, average G's
% average specific energy dissapation per second per square meter 
clc; clear; close all;


% environment constants
mu = 3.5316E12; % gravitational constant (m^3/s^2)
Rb = 600000; % mean radius of body (m)
atmosAlt = 70000; % maximum altitude of atmosphere(m)
planet = 'k'; % determines which planet atmosphere to call on

% input orbital elements
% periapsis is the predicted dip point 
Rp = 630000; % periapsis (m)
Ra = 12000000; % apoapsis (m)

% input craft dimensions
mRocket = 840+500+1000+130+40+300+250+100 % mass of the rocket (kg)
diameter = 1.25; % diameter of the ship (m)
length = 8; % length of the ship (m)
% orientation is the entry angle from the ship longitudinal axis 
% 90 degrees = sideways, 0 degrees = head on
aoa = 0*pi()/180; % angle of attack (radians)
cd = 1; % coefficient of drag (frontal area, see https://www.grc.nasa.gov/www/k-12/airplane/shaped.html)
cl = 2*pi()*aoa; % lift component for wings (see http://brennen.caltech.edu/fluidbook/externalflows/lift/flatplateairfoil.pdf)

% calculate ship areas 
frontal_area = pi()*diameter^2/4; % m^2
lateral_area = diameter*length; % m^2


% simulation constants 
tstep = 0.1; % time steps 
tmax = 3600; % maximum time 
time = [0:tstep:tmax]; % time iterable

% simulation variables
endStep = 0;
maxSimSteps = size(time,2);
% calculate 10% marker for steps 
marker = floor(maxSimSteps / 40);


% vectors: r = radial out, p = tangent to orbit 
position = zeros(maxSimSteps,2); % position over time (r, p)
velocity = zeros(maxSimSteps,2); % velocity over time (r, p)
acceleration = zeros(maxSimSteps,2); % acceleration vectors (r, p)

% scalars (direction is always known)
weight = zeros(maxSimSteps,1); % radial in
lift = zeros(maxSimSteps,1); % lift vector (perpendicular to velocity vector)
drag = zeros(maxSimSteps,1); % negative velocity vector 
q = zeros(maxSimSteps,1); % dynamic pressure, negative velocity vector 

attitude = zeros(maxSimSteps,1); % aircraft attitude

% set up initial position and velocity (at entry)
% get orbital elements 
% calculate entry point 

% get velocity vector in polar coordinates at the intersect point 
% also get true anomaly (ta)
[velocity(1,:), ta] = get_vel(mu, Rp, Ra, atmosAlt + Rb);

% configure intial starting coordinates
position(1,1) = Rb+atmosAlt;
position(1,2) = -ta;

for step = [1:1:maxSimSteps]
  % drag (always inline with velocity vector)
  if(position(step,1) >= Rb && position(step,1) < Rb + atmosAlt)
    % calculate spacecraft frontal area 
    drag_area = frontal_area * cos(aoa) + lateral_area * sin(aoa); % used to calculate drag 
    lift_area = frontal_area * sin(aoa) + lateral_area * cos(aoa); % used to calculate lift
    % calculate dynamic pressure q (Pa/m^2)
    [q(step), M, a] = q_calc(position(step,1)-Rb,norm(velocity(step,:)),planet);
    drag(step) = -q(step) * drag_area * cd;
    lift(step) = q(step) * lift_area * cl;
  elseif(position(step,1) > atmosAlt + Rb + 1000)
    endStep = step;
    disp('Rocket exited atmosphere.')
    break
  elseif(position(step,1) <= Rb)
    endStep = step;
    disp('Rocket collided with the planet.')
    break 
  endif
  
  % gravitational forces, aligned with radial out
  weight(step) = -mu * mRocket * (position(step,1))^(-2);

  % sum forces in r and p directions
  % attitude is held by thrusters for the given amount of steps 
  % otherwise the rocket is allowed to naturally curve over on its own 
  rForce = lift(step) * cos(attitude(step)) + drag(step) * sin(attitude(step)) + weight(step);
  pForce = lift(step) * sin(attitude(step)) + drag(step) * cos(attitude(step));
  
  % calculate change in velocity 
  % F = ma, a*t = DV
  acceleration(step,1) = rForce / mRocket;
  acceleration(step,2) = pForce / mRocket;
  velocity(step,:) = velocity(step,:) + acceleration(step,:).*tstep;
  
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

fprintf('simulation finished.\n')


% calculate acceleration vector length
normAccel = sqrt(acceleration(1:endStep,1).^2+acceleration(1:endStep,2).^2);

% speed 
speed = sqrt(velocity(1:endStep,1).^2+velocity(1:endStep,2).^2);

% calculate power
power = drag(1:endStep).*speed;

% find peak acceleration 
maxAccel = max(normAccel);
maxGs = maxAccel / 9.81; % peak g's
fprintf('Peak acceleration = %.2f m/s^2, %.2f g\n',maxAccel,maxGs)

% integrate drag across time to get total energy loss across the pass 
totalELost = sum(power.*tstep);
specificELost = totalELost / mRocket;
fprintf('Total energy lost = %.4d J, specific energy lost = %.4d J/kg\n',totalELost,specificELost)

% average power loss 
avgPLost = totalELost / time(endStep);
% average specific power loss (mass)
avgSpecificPLost = specificELost / time(endStep);
fprintf('Average power loss = %.4d W, average specific power loss = %0.4d W/kg\n', avgPLost, avgSpecificPLost)

% peak power
peakPower = min(power);
% peak specific power 
peakSpecificPower = peakPower / mRocket;
fprintf('Peak power = %.4d W, peak specific power = %.4d W/kg\n', peakPower, peakSpecificPower)

% specific power per unit area 
peakAreaPower = peakPower / (2*frontal_area + lateral_area);
avgAreaPower = avgPLost / (2*frontal_area + lateral_area);
fprintf('Peak power loss per area = %0.4d W/m^2, average power loss per area = %0.4d W/m^2\n',peakAreaPower,avgAreaPower)

% calculate change in velocity at periapsis through energy loss 

% calculate orbital elements of the original orbit 
[Vp,Va,Period,eccentricity,semilatus_rectum,specificEnergy,angular_momentum] = elliptical_param(mu, Rp, Ra);
newVp = sqrt(2*((Vp^2/2)+specificELost)); % new periapsis velocity
deltaV = newVp - Vp;
fprintf('Delta V at periapsis: %.2f m/s\n',deltaV)
 
% output old orbital parameters  
fprintf('\nOriginal Orbital Parameters: Vp = %.1f m/s, Rp = %.5d m, Va = %.1f m/s, Ra = %0.5d m\n',Vp, Rp, Va, Ra)
fprintf('T = %.0f, e = %.3f, p = %.5d\n',Period,eccentricity,semilatus_rectum)


% calculate new orbital parameters
if(position(endStep,1) > Rb)
[Rpnew,Ranew,Vpnew,Vanew,Tnew,enew,pnew] = get_elliptical_from_state(mu, position(endStep,1), speed(endStep), position(endStep,1) * velocity(endStep,2));
fprintf('New Orbital Parameters: Vp = %.1f m/s, Rp = %.5d m, Va = %.1f m/s, Ra = %0.5d m\n',Vpnew, Rpnew, Vanew, Ranew)
fprintf('T = %.0f, e = %.3f, p = %.5d\n',Tnew,enew,pnew)
else 
fprintf('Full planetary capture: rocket collided with planet going %.1f m/s\n',speed(endStep))
endif 

% Pass statistics 
fprintf('Total entry time = %.1f s\n',time(endStep))
fprintf('Entry true anomaly = %.1f deg, ending true anomaly = %.1f deg, change in true anomaly = %.1f deg\n',...
         position(1,2)*180/pi(),...
         position(endStep,2)*180/pi(),...
         (position(endStep,2)-position(1,2))*180/pi())

% plot the data to the charts

% specific power per unit mass 
%subplot(1,2,1)
%plot(time(1:endStep),power(1:endStep)./mRocket)
%xlabel('time (s)')
%ylabel('specific power (W/kg)') 
%title('specific power loss vs. time')

% specific power per unit area
%subplot(1,2,2)
figure
plot(time(1:endStep),power(1:endStep)./(2*frontal_area + lateral_area))
xlabel('time (s)')
ylabel('specific area power (W/m^2)') 
title('specific area power loss vs. time')

% plot position vs. angle
%subplot(1,2,1)
figure
polar(position(1:endStep,2),position(1:endStep,1))
hold on
% plot planet
polar(linspace(-ta, ta, 100), Rb.*ones(100,1))
% plot planet atmosphere
polar(linspace(-ta, ta, 100), (Rb+atmosAlt).*ones(100,1))
% plot expected periapsis 
polar(0,Rp,'*')
title('position vs. time')
legend('position','planet','atmosphere','exp. periapsis')

% plot altitude vs. time
figure 
plot(time(1:endStep),position(1:endStep,1)-Rb)
title('altitude vs. time')
xlabel('time (s)')
ylabel('altitude (m)')

% plot acceleration vs. time 
figure 
% plot total acceleration 
plot(time(1:endStep-1),normAccel(1:endStep-1))
hold on 
% plot gravity contribution
plot(time(1:endStep-1),mu./position(1:endStep-1,1).^2)
title('acceleration vs. time')
xlabel('time (s)')
ylabel('acceleration (m/s^2)')
legend('total','gravity')

% plot speed vs. time 
figure 
plot(time(1:endStep),speed)
title('speed vs. time')
xlabel('time (s)')
ylabel('speed (m/s)')

