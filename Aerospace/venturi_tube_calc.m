# Venturi tube calculator
# calculates area ratio depending on flow constants
# Damon Printz
# 8/11/2020

clc; clear;
P1 = 20 * 6894.76;  % incoming pressure (Pa)
P2 = 10 * 6894.76; % pressure in the neck (Pa)
V2 = 343;     % speed through the venturi neck (m/s)
rho = 1.225;  % air density (kg/m^3)

if( P2/rho-P1/rho+V2^2/2 < 0 )
  fprintf('cannot achieve desired ratio\n')
else
  % ratio of incoming area A1 to neck area A2  
  A1_A2 = V2*(2*(P2/rho-P1/rho+V2^2/2))^(-1/2)
  % ratio of incoming diameter d1 to venturi neck diameter d2
  d1_d2 = sqrt(A1_A2)
endif 


