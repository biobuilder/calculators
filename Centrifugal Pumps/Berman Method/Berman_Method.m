# Berman Method Centrifugal Pump Design
# follows Design and Performance Analysis of Centrifugal Pump
# by Khin Cho Thin, Mya Mya Khaing, and Khin Maung Aye
# https://pdfs.semanticscholar.org/e351/0bae5253c090888ac43368b1a136c4b1d2b4.pdf 
# second reference 
# https://pdfs.semanticscholar.org/f3fa/8d66fd84b93d968019daf4578c3ba7b1824e.pdf

# Damon Printz
# 7-30-2020

clc; clear; close all;
% input parameters 
Q = 0.015;        % flow rate(m^3/s)
H = 20;           % head (m)
n = 1800;         % rpm

% material constants
rho = 1000; % liquid density (kg/m^3)
g = 9.81;   % gravitational constant (m/s^2)


% - - - - - - - - - - - - - - calculations - - - - - - - - - - - - - - - 

% calculate specific speed 
ns = n*sqrt(Q)/(H^(3/4));  % specific speed
fprintf('specific speed = %.1f\n',ns)







%{
% input parameters 
Q = 1E-6 * 2930;  % flow rate from ml/s to (m^3/s)
H = 10;           % head (m)
n = 2900;         % rpm

% material constants
rho = 1000; % liquid density (kg/m^3)
g = 9.81;   % gravitational constant (m/s^2)


% calculations

ns = 3.65*n*sqrt(Q)/(H^(3/4));  % specific speed
fprintf('specific speed = %.1f\n',ns)

N = rho*g*H*Q;  % water power 
fprintf('water power = %.1f W\n', N)

a1 = 1; % safety factor in charge condition of working pump
n0 = 0.5; % pump efficiency
Mmax = a1 * rho * g * H * Q / n0; % shaft power
fprintf('maximum shaft power = %.1f W\n', Mmax)

% calculate inlet diameter of impeller
K0 = 4.5; % constant chosen for equation
D1 = 1.35 * K0 * sqrt(Q / n);
fprintf('inlet diameter of impeller = %.1f mm\n',D1*1000)

% calculate eye diameter of impeller 
D0 = K0*(Q/n)^(1/3);
fprintf('eye diameter of impeller = %.1f mm\n',D0*1000)

% calculate outlet diameter of impeller 
D2 = 19.2*(ns/100)^(2/6)*sqrt(2*g*H)/n;
fprintf('outlet diameter of impeller = %.1f mm\n',D2*1000)

% calculate shaft diameter at hub section 
T = 9.65*N/n; % torsional moment estimate
t = 75000;  % not sure what this is 
dsh = sqrt(T/(0.2*t));
fprintf('shaft diameter at hub section = %.1f mm\n',dsh*1000)

% hub diameter 
Dbt = (1.75) * dsh;
fprintf('hub diameter = %.1f mm\n', Dbt*1000)

% hub length 
Lbt = 2*dsh;
fprintf('hub length = %.1f mm\n', Lbt*1000)

% inlet width of impeller 
R0 = D0/2; % radius of impeller eye 
b1 = R0/2; % inlet width of impeller 
fprintf('inlet width of impeller = %.1f mm\n', b1*1000)

% outlet width of impeller 
b2 = 0.78*sqrt(ns/100)*(Q/n)^(1/3);
fprintf('outlet width of impeller = %.1f mm\n', b2*1000)

% hydraulic efficiency 
Er = 1-0.42/(log(D0)-0.172)^2;
fprintf('hydraulic efficiency = %f\n', Er)

% leakage head 
HT = H/Er; % pressure head 
Dy = D0+10; % seal diameter 
Hy = HT-V2^2/(2*g)-U2^2/(8*g)*(1-Dy^2/D2^2);
fprintf('leakage head = %f\n',Hy)

% minimum clearance between wear ring and casing
del = Dy/1000;
fprintf('min clearance between wear ring and casing = %f mm\n', del*1000)

% inlet blade angle of impeller 
beta1 = atan(U1/Vm1);
fprintf('inlet blade angle of impeller = %f deg\n',beta1*180/pi())

% outlet blade angle of impeller 
beta2 = 40*pi()/180;
fprintf('outlet blade angle of impeller = %.1f deg\n',beta2*180/pi())

% blade number 
Z = 6.5*(D2+D1)/(D2-D1)*sin(beta1+beta2)/2;
fprintf('blade number = %f\n',Z)
%}









