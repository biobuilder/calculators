% Attempts to find a guassian trajectory between two vectors 
% with a specified transit time 
% Damon Printz
% 4/5/2020

function [ hyperExcess1, ...      % departure hyperbolic excess velocity 
           hyperExcess2, ...      % arrival hyperbolic excess velocity 
           inc, ...               % trajectory inclination [rad]
           aop, ...               % trajectory argument of periapsis [rad]
           lan, ...               % trajectory longitude of ascending node [rad]
           e, ...                 % trajectory eccentricity
           lpt ...                % launch point ascending node [rad]
           ] = solve_gaussian_trajectory (...
           r1, ...                % departure vector 
           r2, ...                % arrival vector 
           mu, ...                % gravitational parameter
           tt, ...                % transit time 
           Vp1, ...               % departure planet velocity 
           Vp2, ...               % arrival planet velocity 
           timeAccu, ...          % trajectory accuracy 
           maxTimeIter, ...       % maximum iterations
           maxDV ...              % maximum Delta V allowed 
           )  
  % debugging flag 
  verbose = 0;  
  
  % predefine all outputs  
  hyperExcess1 = 0;
  hyperExcess2 = 0;
  inc = 0;
  aop = 0;
  lan = 0;
  e = 0;
  lpt = 0;
    
  % calculate the length of the departure and arrival vectors     
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

  % initialize guess variables for determining transfer orbital parameters 
  pnew = 0;
  p = 0;
  t = 0;
  po = 0;
  to = 0;

  trajCnt = 1;
  while( abs(t - tt) > timeAccu && trajCnt < maxTimeIter)
    
    % select trial value for semilatus rectum
    if( trajCnt == 1 )
      p = pstart1;
    elseif( trajCnt == 2 )
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
    
    % make sure division is usable 
    if(a*(1-f) == 0)
      if(verbose != 0)
        fprintf('sg: division by zero\n')
      endif
      return
    endif 
    
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
    trajCnt = trajCnt + 1;
  endwhile

  % output calculated values 
  if(verbose != 0)
    fprintf('sg: p = %.0f m, t = %.0f s',p,t)
  endif 
  if(trajCnt < maxTimeIter)
    if(verbose != 0)
      fprintf(' in %i iterations\n',trajCnt)
    endif 
  else 
    if(verbose != 0)
      fprintf('sg: - did not converge\n')
    endif 
    return
  endif 

  % error checking 
  if(Dtanom < pi() && p > p1 && trajCnt < maxTimeIter)
    if(verbose != 0)
      fprintf('sg: valid elliptical solution, p > p1\n')
    endif 
  elseif(Dtanom > pi() && p < p2 && trajCnt < maxTimeIter)
    if(verbose != 0)
      fprintf('sg: valid hyperbolic solution, p < p2\n')
    endif 
  else 
    % return because solution was not valid 
    if(verbose != 0)
      fprintf('sg: invalid solution\n')
    endif
    return
  endif 

  
  
  % Calculate departure and intercept velocity vectors - - - - - - - - - - - - - -
  if(p*(1-cos(Dtanom)) == 0)
    if(verbose != 0)
      fprintf('sg: r1 == r2 already\n')
    endif 
    return
  endif 
  
  gdot = 1-r1len/p*(1-cos(Dtanom));
  
  v1 = (r2 - f.*r1)/g;
  v2 = fdot.*r1+gdot.*v1;

  if(abs(v1) + abs(v2) > maxDV)
    if(verbose != 0)
      fprintf('sg: max DV exceeded\n')
    endif 
    return
  endif 
  
  % print off velocities 
  if(verbose != 0)
    fprintf('sg: v1 = %.1f m/s, v2 = %.1f m/s\n',norm(v1),norm(v2))
  endif 
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
  inc = acos(hvect(3) / h);

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
  lpt = acos(dot(nvect,r1)/(n*r1len));

  % display orbital elements 
  if(verbose != 0)
    fprintf('sg: orbital elements: e = %.5f, i = %.3f deg, lan = %.1f deg, aop = %.1f deg\n', ...
            e, inc*180/pi(), lan*180/pi(), aop*180/pi())
    % display launch information 
    fprintf('sg: initial true anomaly = %.3f deg, launch point (ascending node) %.1f deg\n',...
            tanom1*180/pi(), lpt*180/pi())
  endif 


  % Calculate hyperbolic excess velocity at departure, injection DV, - - - - - - -
  % and zenith angle of departure 

  % Hyperbolic Departure - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  % 1. Escape from host body sphere of influence

  % calculate excess hyperbolic escape velocity from starting point 
  hyperExcessVect1 = v1 - Vp1;
  hyperExcess1 = norm(hyperExcessVect1);
  
  hyperExcessVect2 = v2 - Vp2;
  hyperExcess2 = norm(hyperExcessVect2);

  if(verbose != 0)
    fprintf('sg: departure hyperbolic excess escape velocity = %.1f m/s\n',hyperExcess1)
    fprintf('sg: arrival hyperbolic excess escape velocity = %.1f m/s\n',hyperExcess2)
  endif 
  % calculate injection velocity 
  %injV1 = sqrt(hyperExcess1^2+2*mu_host/r0len);
  %fprintf('orbital injection velocity from surface = %.1f m/s\n', injV1)

  % calculate injection DV to orbit 
  %DV1 = injV1 - sqrt(mu_host / r0len);
  %fprintf('departure DV from orbit = %.1f m/s\n', DV1)

  % calculate zenith angle of departure 
  %zenithAOD = acos( dot(r1,hyperExcessVect1) / (r1len*norm(hyperExcessVect1)) );
  %fprintf('zenith angle of departure = %.3f deg\n', zenithAOD*180/pi())

  % get orbital elements of the hyperbolic orbit
  %ahyp1 = -mu_host/hyperExcess1^2; % semimajor axis 
  %ehyp1 =  -r0len/ahyp1+1; % eccentricity 
  %phyp1 = ahyp1*(ehyp1^2-1); % semilatus rectum 
  %bhyp1 = -ahyp1 * sqrt(ehyp1^2-1); % impact parameter (semiminor axis)

endfunction
