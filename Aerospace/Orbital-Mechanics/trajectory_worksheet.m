% Trajectory Analysis Worksheet
% Reference: http://www.braeunig.us/space/interpl.htm
% Damon Printz
% 4/3/2020

clc; clear; close all;

% Input simulation parameters - - - - - - - - - - - - - - - - - - - - - - - - -

% gravitational parameter of the sun Kerbol
mu = 1.327124E20; %*(149.59787E9)^(-3); %1.1723328E18; 

% note: r1 and r2 cannot be colinear
% starting vector [m]
r1 = [.473265, -0.899215, 0].*(149.59787E9)^(1); %[13599840256, 0 , 0];

% ending vector [m]
r2 = [0.066842, 1.561256, 0.030948].*(149.59787E9)^(1); %[9931011387*cosd(175), 9931011387*sind(175), 0];

% transit time (s)
tt = 207*24*3600; %2 * pi() * sqrt( (a^3) / mu ) * 0.5; % 5657995 * (3/4);

TIME_ACCURACY = 3600; % must be within 1 hour of target time
MAX_TIME_ITERATIONS = 100; % maximum tries to get the transit time accurate 



% target planet sphere of influence [m]
%SOI = 85109365;

% host body gravitational parameter 
mu_host = 3.986005E14;
% host body velocity vector at departure r1
Vp = [25876.6, 13759.5, 0]; 
% circular orbital radius around host
r0len = (6378.14+200)*1000;


% target body gravitational parameter 
mu_target = 4.282837E13;
% target body velocity vector at approach r2
% circular orbital velocity around target 





% Search for valid trajectory - - - - - - - - - - - - - - - - - - - - - - - - -

% calculate length of vectors r1, r2
r1len = norm(r1);
r2len = norm(r2);

% calculate change in true anomaly between the two vectors
Dtanom = acos(dot(r1,r2)/(r1len*r2len));

% set intermediate constants
k = r1len * r2len * (1-cos(Dtanom));
l = r1len + r2len;
m = r1len*r2len*(1+cos(Dtanom));

p1 = k / (l + sqrt(2*m));
p2 = k / (l - sqrt(2*m));

% check to see if change in true anomaly is short way (less than pi())
if(Dtanom < pi())
  % short way
  % we know p1 < p < inf
  pstart1 = p1*1.1;
  pstart2 = p1*1.5;
else 
  % long way 
  % we know 0 < p < p2
  pstart1 = p2 * 0.6;
  pstart2 = p2 * 0.4;
endif 
  
pnew = 0;
p = 0;
t = 0;
po = 0;
to = 0;

i = 1;
while( abs(t - tt) > TIME_ACCURACY && i < MAX_TIME_ITERATIONS)
  
  % select trial value for semilatus rectum
  if(i == 1)
    p = pstart1;
  elseif( i == 2 )
    p = pstart2;
  else
    p = pnew;
  endif 

  % calculate a parameter
  a = m*k*p / ((2*m-l^2)*p^2+2*k*l*p-k^2);
  
  % calculate f, g, fdot
  f = 1-r2len/p*(1-cos(Dtanom));
  g = r1len*r2len*sin(Dtanom)/sqrt(mu*p);
  fdot = sqrt(mu/p)*tan(Dtanom/2)*((1-cos(Dtanom))/p-1/r1len-1/r2len);
  
  % elliptical orbits have a > 0
  if(a > 0) 
    DE = acos(1-r1len/a*(1-f));
    t = g+sqrt(a^3/mu)*(DE-sin(DE));
  else
    % hyperbolic orbits have a < 0
    % t is calculated differently 
    DF = acosh(1-r1len/a*(1-f));
    t = g+sqrt(abs(a)^3/mu)*(sinh(DF)-DF);
  endif 
  
  % calculate new trial value for semilatus rectum
  pnew = p + (tt - t) * (p - po) / (t - to);
  
  % record the old values 
  po = p;
  to = t;
  i = i + 1;
endwhile

% output calculated values 
fprintf('p = %.0f m, t = %.0f s',p,t)
if(i < MAX_TIME_ITERATIONS)
  fprintf(' in %i iterations\n',i)
else 
  fprintf(' - did not converge\n')
endif 

% error checking 
if(Dtanom < pi() && p > p1)
  fprintf('valid elliptical solution, p > p1\n')
elseif(Dtanom > pi() && p < p2)
  fprintf('valid hyperbolic solution, p < p2\n')
else 
  p = 0;
  fprintf('invalid solution\n')
endif 


if(p > 0)

% Calculate departure and intercept velocity vectors - - - - - - - - - - - - - -
gdot = 1-r1len/p*(1-cos(Dtanom));

v1 = (r2 - f.*r1)/g;
v2 = fdot.*r1+gdot.*v1;

% print off velocities 
fprintf('v1 = %.1f m/s, v2 = %.1f m/s\n',norm(v1),norm(v2))

% Calculate orbital elements - - - - - - - - - - - - - - - - - - - - - - - - - -
hvect = cross(r1,v1); % get specific angular momentum vector 
h = norm(hvect); % get total specific angular momentum 

% "line of nodes" n is perpendicular to Z axis and h vector 
% this node points in the direction of the ascending node 
nvect = [-hvect(2),hvect(1),0]; % get line of nodes 
n = norm(nvect);

% calculate direction of periapsis
% the magnitude of which is the eccentricity of the orbit
evect = ((norm(v1)^2-mu/r1len).*r1 - dot(r1,v1) .* v1 ) / mu;
e = norm(evect);

% calculate inclination
i = acos(hvect(3) / h);

% calculate longitude of ascending node 
lan = 0;
if(nvect(2) > 0)
  lan = acos(nvect(1)/n);
else 
  lan = 2*pi() - acos(nvect(1)/n);
endif 

% calculate argument of periapsis
aop = 0;
if(evect(3) > 0)
  aop = acos(dot(nvect,evect)/(n*e));
else 
  aop = 2*pi()-acos(dot(nvect,evect)/(n*e));
endif 

% initial true anomaly 
tanom1=0;
if(dot(r1,v1) > 0)
  tanom1 = acos(dot(evect,r1)/(e*r1len));
else 
  tanom1 = 2*pi() - acos(dot(evect,r1)/(e*r1len));
endif 

% launch point (ascending node) 
asnode = acos(dot(nvect,r1)/(n*r1len));

% display orbital elements 
fprintf('orbital elements: e = %.5f, i = %.3f deg, lan = %.1f deg, aop = %.1f deg\n', ...
        e, i*180/pi(), lan*180/pi(), aop*180/pi())
% display launch information 
fprintf('initial true anomaly = %.3f deg, launch point (ascending node) %.1f deg\n',...
        tanom1*180/pi(), asnode*180/pi())


% Calculate hyperbolic excess velocity at departure, injection DV, - - - - - - -
% and zenith angle of departure 

% Hyperbolic Departure - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% 1. Escape from host body sphere of influence

% calculate excess hyperbolic escape velocity from starting point 
hyperExcessVect1 = v1 - Vp;
hyperExcess1 = norm(hyperExcessVect1);
fprintf('departure hyperbolic excess escape velocity = %.1f m/s\n',hyperExcess1)

% calculate injection velocity 
injV1 = sqrt(hyperExcess1^2+2*mu_host/r0len);
fprintf('orbital injection velocity from surface = %.1f m/s\n', injV1)

% calculate injection DV to orbit 
DV1 = injV1 - sqrt(mu_host / r0len);
fprintf('departure DV from orbit = %.1f m/s\n', DV1)

% calculate zenith angle of departure 
zenithAOD = acos( dot(r1,hyperExcessVect1) / (r1len*norm(hyperExcessVect1)) );
fprintf('zenith angle of departure = %.3f deg\n', zenithAOD*180/pi())

% get orbital elements of the hyperbolic orbit
ahyp1 = -mu_host/hyperExcess1^2; % semimajor axis 
ehyp1 =  -r0len/ahyp1+1; % eccentricity 
phyp1 = ahyp1*(ehyp1^2-1); % semilatus rectum 
bhyp1 = -ahyp1 * sqrt(ehyp1^2-1); % impact parameter (semiminor axis)

% 2. Enter parking orbit around target



endif
