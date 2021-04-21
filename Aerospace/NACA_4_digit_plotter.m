% plots NACA 4-digit airfoils and prints the coordinates to the command window
% based on Wikipedia Article https://en.wikipedia.org/wiki/NACA_airfoil
% Damon Printz
% 4/11/2021

clc; clear; close all;

% resolution of the plotter
points = 15;

% maximum thickness as percentage of the cord (digits 3 and 4)
t = 0.60;

% location of maximum camber as percentage of chord (digit 2)
p = 0.30;

% maximum camber (digit 1)
m = 0; 


% set up equally spaced coordinates along the chord from 0 to 1
x = linspace(0,1,points);

% calculate half chord width
yt = (5*t).*(0.2969.*sqrt(x)-0.1260.*x-0.3516.*(x.^2)+0.2843.*(x.^3)-0.1015.*(x.^4));

yc = zeros(size(yt));
ycCnt = 1;
% calculate mean camber line 
for xVal = x(x<p)
  yc(ycCnt) = m/p^2*(2*p*xVal - xVal^2);
  ycCnt = ycCnt + 1;
endfor
for xVal = x(x>=p)
  yc(ycCnt) = m/(1-p)^2*((1-2*p)+2*p*xVal-xVal^2);
  ycCnt = ycCnt + 1;
endfor

% apply perpendicular to camber line 
dycdx = zeros(size(yt));
dycCnt = 1;
for xVal = x(x<p)
  dycdx(dycCnt) = 2*m/p^2*(p-xVal);
  dycCnt = dycCnt + 1;
endfor
for xVal = x(x>=p)
  dycdx(dycCnt) = 2*m/(1-p^2)*(p-xVal);
  dycCnt = dycCnt + 1;
endfor

theta = atan(dycdx);


% upper coordinate list 
xu = x - yt .* sin(theta);
yu = yc + yt .* cos(theta);

% lower cooridinate list 
xl = x + yt .* sin(theta);
yl = yc - yt .* cos(theta);

% print the coordinates 

fprintf('upper coordinates, (x, y) pairs\n')
for ind = [1:1:points]
  fprintf('%f,%f\n',xu(ind),yu(ind))
endfor

fprintf('lower coordinates, (x, y) pairs\n')
for ind = [1:1:points]
  fprintf('%f,%f\n',xl(ind),yl(ind))
endfor



% plot the airfoil
plot(xu,yu)
hold on 
plot(xl,yl)
axis("equal")
grid on 