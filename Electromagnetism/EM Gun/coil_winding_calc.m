% EM Gun coil calculations
% calculates coil winding characteristics
% Damon Printz
% 12/26/2020

clc; close all; clear;

%{
fprintf('Wound copper wire stack\n')
wire_dia = 0.511 % mm, wire diameter

coil_id = 12 % mm, inner diameter of coil

layers = [1:1:10]; % number of layers to analyze 

% material constants
wire_resistivity = 1.724E-8; % ohm m, copper
permeability = 4*pi()*10^(-7); % T/amp m, magnetic permeability of free space
k = 1; % relative permeability of inner core material (air)
%k = 200; % relative permability of iron
wire_density = 8.96; % g/cm^3, copper 

% calculations
wraps_per_m = 1000/wire_dia % number of wraps per m per layer
wire_area = (pi()*wire_dia^2)/4; % mm^2, cross-sectional area
wire_resistance_per_m = wire_resistivity * 1000000 / wire_area % ohms per m

% magnetic field strength per layer per amp
B = permeability * k * wraps_per_m;
fprintf('B = %f T / amp per layer\n',B)

% wire length for the given range of layers 
fprintf('layers, wire length per meter of barrel, resistance, mass, miliTeslas/amp, outer diameter\n')
wireLen = 0;
wireMass = 0;
for layer = layers
  wireLen = wireLen + wraps_per_m * (coil_id+(layer-0.5)*wire_dia) * pi() / 1000; % length in m
  wireMass = (wireLen * wire_area) * wire_density / 1000; % mass in kg
  wireResistance = wireLen * wire_resistance_per_m; % resistance in Ohms
  magStrength = B*layer*1000; % in miliTeslas 
  coilOD = coil_id + layer * wire_dia; % outer diameter of coil
  fprintf('%i, %.0f m, %.1f Ohm, %.2f kg, %.2f mT/A, %.0f mm\n', layer, wireLen, wireResistance, wireMass, magStrength, coilOD) 
endfor
%}


%
fprintf('PCB wire stack\n')
wire_wid = 0.1524 % mm, wire width
wire_layer_wid = 0.3048; % mm, width of each layer radially out  
wire_thick = 0.03556 % mm, wire thickness
coil_id = 12 % mm, inner diameter of coil
wraps_per_m = 4/0.4*1000 % number of wraps per m per layer

coils_per_pcb = 4 % number of connected coils per pcb 

layers = [10:1:100]; % number of layers to analyze 

% material constants
wire_resistivity = 1.724E-8; % ohm m, copper
permeability = 4*pi()*10^(-7); % T/amp m, magnetic permeability of free space
k = 1; % relative permeability of inner core material (air)
%k = 200; % relative permability of iron
wire_density = 8.96; % g/cm^3, copper 

% calculations
wire_area = wire_wid * wire_thick; % mm^2, cross-sectional area
wire_resistance_per_m = wire_resistivity * 1000000 / wire_area % ohms per m

% magnetic field strength per layer per amp
B = permeability * k * wraps_per_m;
fprintf('B = %f T / amp per layer\n',B)

% wire length for the given range of layers 
fprintf('layers, wire length per meter of barrel, resistance, mass, miliTeslas/amp, outer diameter\n')
wireLen = 0;
wireMass = 0;
for layer = layers
  wireLen = wireLen + wraps_per_m * (coil_id+(layer-0.5)*wire_layer_wid) * pi() / 1000; % length in m
  wireMass = (wireLen * wire_area) * wire_density / 1000; % mass in kg
  %wireResistance = wireLen * wire_resistance_per_m; % resistance in Ohms
  % PCB coils may need to be connected in series to increase resistance 
  % coils are connected in pairs or quadruplets 
  wireResistance =  (coils_per_pcb * (wireLen / wraps_per_m) * wire_resistance_per_m)/wraps_per_m; % ohms per single coil
  magStrength = B*layer*1000; % in miliTeslas 
  coilOD = coil_id + layer * wire_layer_wid; % outer diameter of coil 
  if(mod(layer, 10) == 0)
    fprintf('%i, %.0f m, %.5f Ohm, %.2f kg, %.2f mT/A, %.0f mm\n', layer, wireLen, wireResistance, wireMass, magStrength, coilOD) 
  endif
endfor
%
