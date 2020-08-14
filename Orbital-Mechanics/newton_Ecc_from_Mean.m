% Newton's Method solver for getting Eccentric anomaly from Mean anomaly
% Damon Printz
% 4/5/2020


% Inputs
% M = mean anomaly [rad]
% e = eccentricity
% minimum accuracy = accu (angle change) [rad]
% maxIter = maximum iterations 

% Outputs
% E = eccentric anomaly [rad] 

function E = newton_Ecc_from_Mean (M, e, accu, maxIter)
% minimum angle change accu 
% accu = 0.01*pi()/180; % [rad]
%maxIter = 20; % maximum iterations 
startVal = M;

cnt = 0; % current iteration count 
delta = accu+1; % angle change 

xold = startVal;
while(cnt < maxIter && abs(delta) > accu)
  
  xnew = xold - (M - (xold + e*sin(xold))) / (-(1+e*cos(xold)));
  
  delta = xnew - xold;
  xold = xnew;
  cnt = cnt + 1;
endwhile 

E = xnew;

%fprintf('E = %.5f after %i iterations\n',E,cnt)

if(delta > accu)
fprintf('Warning: angle did not converge\n')
endif 

endfunction
