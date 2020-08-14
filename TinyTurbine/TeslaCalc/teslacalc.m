% Tiny Turbine - Tesla Turbine Calculations
% Damon Printz
% 2/17/2018

% additional references:
% http://homepages.engineering.auckland.ac.nz/~pkel015/SolidMechanicsBooks/Part_II/04_ElasticityPolar/ElasticityPolars_04_BodyForcesRotatingDiscs.pdf
% simplified calculations for different situations
% http://www.roymech.co.uk/Useful_Tables/Mechanics/Rotating_cylinders.html

clear; clc; close all;
more off % keep command prompt from prompting too much

% disk dimensions
R = 6.0E-3; % outer radius (m)
R0 = 0.0E-3; % inner radius (m)

% disk rpm
rpm = 10000; % revolutions per minute
w = 2*pi*rpm/60; % radians per second
% tip velocity
tip_vel = w*R;

% select material data file to load
% 'AISI_1020_mild_steel.m'
% 'steel_1144.m'
% 'aluminum_6061_T6.m'
% 'brass_UNS_C36000.m'
% 'abs_generic.m'
% 'acrylic_generic.m'
folder = 'mat_prop/';
material = 'brass_UNS_C36000.m';
mat_data = load(strcat(folder,material));
rho = mat_data(1);
v = mat_data(2);
E = mat_data(3);
fprintf('%s properties:\n',material(1:length(material)-2))
fprintf('  yield strength = %f MPa\n', mat_data(4)*1E-6)
fprintf('  shear strength = %f MPa\n', mat_data(6)*1E-6)

% General solution for a rotating ring
% http://www.amesweb.info/StructuralAnalysisBeams/Stresses-Rotating-Rings.aspx
% assumptions: effect of deformation can be ignored
% ring is rigid

r = linspace(R0, R, 100); % create an array of radii for calculations
% radial stress as a function of radius
sigma_r = (3+v)/8*rho*w^2.*(R^2+R0^2-R^2*R0^2./r.^2-r.^2);
% tangental stress as a function of radius
sigma_t = rho*w^2/8*((3+v).*(R^2+R0^2+R^2*R0^2./r.^2)-(1+3*v).*r.^2);

% maximum radial stress
max_rad_stress = (3+v)/8*rho*w^2*(R-R0)^2;
% max tangental stress
max_tan_stress = rho*w^2/4*((3+v)*R^2+(1-v)*R0^2);

% dimensional change in Outer Radius
delta_R = rho*w^2*R/(4*E)*((1-v)*R^2+(3+v)*R0^2);
% dimensional change in Inner Radius
delta_R0 = rho*w^2*R0/(4*E)*((3+v)*R^2+(1-v)*R0^2);

fprintf('max radial stress = %.3f MPa\n',max_rad_stress*1E-6)
fprintf('max tangental stress = %.3f MPa\n',max_tan_stress*1E-6)
fprintf('outer radius change = %.2f mm\n',delta_R*1000)
fprintf('inner radius change = %.2f mm\n',delta_R0*1000)
fprintf('tip velocity = %f m/s\n',tip_vel)

plot(r.*1000, sigma_r.*1E-6)
hold on
plot(r.*1000, sigma_t.*1E-6)
title('stress vs. radial location')
xlabel('radius(mm)')
ylabel('stress (MPa)')
legend('radial','tangental')