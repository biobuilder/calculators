% Hydrofoil Engine Sizing
% David Larson and Damon Printz
% 5/9/2018

clear; clc; close all;
% unit conversion factors
kw_to_hp = 0.7457;    % kw to hp
lb_to_kg = 0.453592;  % lbs to kg
mph_to_mps = 0.44704; % mph to m/s
lb_to_kg = 0.453592;  % lbs to kg


% flight conditions
speed = 30;             % speed in mph
v = mph_to_mps*speed;  % velocity in m/s

% craft specifications
m_pilot = 200;      % pilot mass (lbs)
m_frame = 30;       % frame mass (lbs)
m_engine = 60;      % engine mass (lbs)
m_total = (m_pilot + m_frame + m_engine)*lb_to_kg; % total mass (kg)

b = 255*2 / 1000;   % span in m
S = b * 155 / 1000; % surface area in m^2
AR = b^2 / S;       % aspect ratio

e = 0.8;            % span efficiency factor

% engine specifications
n_engine = 0.5;     % efficiency of system (nominal to mechanical power)

% nature constants
g = 9.81;           % gravitational constant (m/s^2)
rho_w = 997;        % density of pure water (kg/m^3)

% calculate lift and drag - - - - - - - - - - - - - - - - - - - - - - - - - - -
q = 0.5 * rho_w * v^2;    % dynamic pressure (Pa) 
L = g * m_total;          % required lift (N)
%Di = L^2/(q*pi*b^2);     % induced drag, planar wing (N)


CL = L/(q*S);             % coefficient of lift
CDi = CL^2/(pi*e*AR);     % coefficient of induced drag

Di = q*S*CDi;             % induced drag

% calculate flow power
Preq_i = Di * v / 1000;   % required power in kw (from induced drag)

% required engine power (nominal, kw)
Pnom = Preq_i / n_engine; 

% display power requirements
fprintf('Theoretical power requirement = %.2f hp\n',Preq_i*kw_to_hp)
fprintf('Nominal power requirement: %.2f Hp\n',Pnom*kw_to_hp)
fprintf('Engine efficiency = %.2f\n',n_engine)
