%   Electromagnetic Coil Power Calculator
%   Calculates the power produced by a coil moving through a changing magnetic field
%   Coil coordinates are from the center of the coil
%   Damon Printz
%   7/3/2018

clc; clear; close all;
% coil dimensions (rectangular type)
coil_N = 99; % number of turns of wire
wire_awg = 18; % wire gauge

coil_W = 1; % innner width of the coil, perpendicular to motion [mm]
coil_L = 15; % inner length of the coil, parallel with motion [mm]

coil_H = 10; % maximum coil height (mm)
% coil material properties
wire_resistivity = 1.7E-8; % ohms/meter
% wire outputs
wire_d = 0.127 * 92^((36-wire_awg)/39); % calculate wire diameter from gauge (mm)

n_wraps_axial = floor(coil_H/wire_d); % how many wraps moving upward axially
n_wraps_radial = ceil(coil_N/n_wraps_axial); % how many wraps moving outward radially

coil_T = n_wraps_radial*wire_d*sin(pi/3); % hexagonally packed thickness of coil

coil_T_to_H = coil_T / coil_H; % thickness to height ratio
fprintf('coil wrap thickness/height = %.2f\n',coil_T_to_H)

wire_length = 2*coil_N*(coil_L+coil_T+coil_W); % get coil wire length
coil_R = wire_resistivity*1000*wire_length*4/(pi*wire_d^2);
fprintf('wire length = %.1f m\n',wire_length/1000)
fprintf('coil resistance = %f ohms\n',coil_R)

% The magnetic field density can be described by a sinusoid when magnets
% of alternating polarity are stacked one after the other. The coil can
% be positioned between these magnets to take advantage of the magnetic
% field.  Below are the magnetic field function constants.
B_amplitude = 0.3;   % maximum flux density (T)
B_offset = 0;        % position offset (mm)
magnet_length = 20;  % single magnet length plus spacing to the next one (mm)
% The magnetic flux through the coil rectangle centered at coordinate x
% is equivalent to the following function:
% phi_B = -coil_W*B_amplitude*magnet_length/pi*cos(pi/magnet_length*x){from x - coil_L/2 to x + coil_L/2)

% If x(t) is a sinusoidal function, given by these characteristics:
x_freq = 1;      % frequency of the response (Hz)
x_amplitude = 10; % amplitude of motion (mm)
x_offset = 0;     % offset from original equation

x_period = 1/x_freq;  % period of the response (sec)
t = [0:x_period/100:x_period*2]; % t values
x_t = x_amplitude.*sin(t.*(2*pi/x_period)); % position vs. time
dx_dt = x_amplitude.*(x_period/(2*pi)).*cos(t.*(2*pi/x_period)); % velocity vs. time
x1_t = x_t + coil_L/2; % rightmost coil edge coordinate
x0_t = x_t - coil_L/2; % leftmost coil edge coordinate

% flux with respect to time
phi_B = -coil_W*B_amplitude*magnet_length/pi.*(cos(pi/magnet_length.*x1_t)-cos(pi/magnet_length.*x0_t));
% flux changing with respect to time
dphiB_dt = (coil_W*B_amplitude).*dx_dt.*(sin((pi/magnet_length).*x1_t)-sin((pi/magnet_length).*x0_t));
EMF = -coil_N .* dphiB_dt;

figure 
subplot(3,1,1)
plot(t,x_t)
title('coil position vs. time')
xlabel('t (sec)')
ylabel('x (mm)')
grid

subplot(3,1,2)
plot(t, phi_B);
title('Flux through coil vs. time')
xlabel('t (sec)')
ylabel('Flux (uWb)')
grid

subplot(3,1,3)
plot(t,EMF)
title('Voltage vs. time')
xlabel('t (sec)')
ylabel('Voltage (uV)')
grid


%{

% magnetic field inputs
% Describes the static magnetic field in terms of a flux density per unit width
% Surprisingly, magnets stacked one on top of the other in alternating pairs
% actually mimic a sinusoidal function
B_amplitude = 0.3; % maximum flux density per unit width (T/mm)
B_phase = 0; % position offset (mm)
B_period = 35;  % width between any magnet pair

% The static magnetic field function appears as follows:
% B_func = B_amplitude * cos(2*pi*(x-B_phase)/B_period);

% The net magnetic flux through any surface is a function of its position
% within the repeating magnetic field, and is equal to the closed loop
% integral of the rectangular surface and the field normal to the surface
% The coil "center" coordinate is at 
x_0 = -0.5 * coil_W;
x_1 = 0.5 * coil_W;
B_int = B_amplitude * B_period / (2*pi) * (sin(2*pi*x_1/B_period)-sin(2*pi*x_0/B_period));

magnetic_flux = B_int * coil_W * 0.001; % get flux density in Webers

fprintf('Total magnetic flux = %f [Wb]\n',magnetic_flux)

% now calculate EMF as a function of t
% If position is described by a sine wave, where x(t) = A*sin(t/Period), 
% then dx/dt = A * Period * cos(t/Period)
% 

%}










