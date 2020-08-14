% Kerbal Space Program Space Parameters
% Damon Printz
% 3/3/2020


% Sun (Kerbol)
mu_Sun = 1.1723328E18; % m^3/s^2
Tsr_Sun = 432000; % sidereal day (s)
Req_Sun = 261600000; % equatorial radius (m)

sun = struct("mu", mu_Sun, "Tsr", Tsr_Sun, "Req", Req_Sun);

% Kerbin Constants - - - - - - - - - - - - - - - - - - - -
mu_Ker = 3.5316E12; % m^3/s^2
Tsr_Ker = 21549.425; % sidereal day (s)
Req_Ker = 600000; % equatorial radius (m)

% solar parameters
Ra_Ker = 13599840256; % apoapsis (m)
Rp_Ker = 13599840256; % periapsis (m)
e_Ker = 0; % orbital eccentricity (rad)
inc_Ker = 0; % orbital inclination (rad)
argp_Ker = 0; % argument of periapsis (rad)
lonasc_Ker = 0; % longitude of ascending node (rad)
Torb_Ker = 9203545; % orbital period (s)
MA_Ker = 3.14; % mean anomaly (rad) at 0s UT

% Kerbin Moons

% Mun Constants
mu_Mun = 6.5138398E10; % m^3/s^2
Tsr_Mun = 138984.38; % sidereal day (s)
Req_Mun = 200000; % equatorial radius (m)

% Mun-kerbin parameters
Ra_Mun = 12000000; % apoapsis (m)
Rp_Mun = 12000000; % periapsis (m)
e_Mun = 0; % orbital eccentricity (rad)
inc_Mun = 0; % orbital inclination (rad)
argp_Mun = 0; % argument of periapsis (rad)
lonasc_Mun = 0; % longitude of ascending node (rad)
Torb_Mun = 9203545; % orbital period (s)
MA_Mun = 1.7; % mean anomaly (rad) at 0s UT

% Minmus Constants
mu_Min = 1.7658000E9; % m^3/s^2
Tsr_Min = 40400; % sidereal day (s)
Req_Min = 60000; % equatorial radius (m)

% Minmus-kerbin parameters
Ra_Min = 47000000; % apoapsis (m)
Rp_Min = 47000000; % periapsis (m)
e_Min = 0; % orbital eccentricity (rad)
inc_Min = 6*pi()/180; % orbital inclination (rad)
argp_Min = 38*pi()/180; % argument of periapsis (rad)
lonasc_Min = 78*pi()/180; % longitude of ascending node (rad)
Torb_Min = 1077311; % orbital period (s)
MA_Min = 0.9; % mean anomaly (rad) at 0s UT


% Duna Constants - - - - - - - - - - - - - - - - - - - - - - 
mu_Dun = 3.0136321E11; % m^3/s^2
Tsr_Dun = 65517.859; % sidereal day (s)
Req_Dun = 320000; % equatorial radius (m)

% solar parameters
Ra_Dun = 21783189163; % apoapsis (m)
Rp_Dun = 19669121365; % periapsis (m)
e_Dun = 0.051; % orbital eccentricity (rad)
inc_Dun = 0.06*pi()/180; % orbital inclination (rad)
argp_Dun = 0*pi()/180; % argument of periapsis (rad)
lonasc_Dun = 135.5*pi()/180; % longitude of ascending node (rad)
Torb_Dun = 19645697.3; % orbital period (s)
MA_Dun = 3.14; % mean anomaly (rad) at 0s UT

% Duna Moons 

% Ike
mu_Ike = 1.8568369E10; % m^3/s^2
Tsr_Ike = 65517.862; % sidereal day (s)
Req_Ike = 130000; % equatorial radius (m)

% Ike-Duna parameters
Ra_Ike = 3296000; % apoapsis (m)
Rp_Ike = 3104000; % periapsis (m)
e_Ike = 0.03; % orbital eccentricity (rad)
inc_Ike = 0.2*pi()/180; % orbital inclination (rad)
argp_Ike = 0*pi()/180; % argument of periapsis (rad)
lonasc_Ike = 0*pi()/180; % longitude of ascending node (rad)
Torb_Ike = 65766.7; % orbital period (s)
MA_Ike = 1.7; % mean anomaly (rad) at 0s UT

% Eve Constants - - - - - - - - - - - - - - - - - - - - - - -
mu_Eve = 8.1717302E12; % m^3/s^2
Tsr_Eve = 805000; % sidereal day (s)
Req_Eve = 700000; % equatorial radius (m)

% solar parameters
Ra_Eve = 9931011387; % apoapsis (m)
Rp_Eve = 9734357701; % periapsis (m)
e_Eve = 0.01; % orbital eccentricity (rad)
inc_Eve = 2.1*pi()/180; % orbital inclination (rad)
argp_Eve = 0*pi()/180; % argument of periapsis (rad)
lonasc_Eve = 15*pi()/180; % longitude of ascending node (rad)
Torb_Eve = 5657995; % orbital period (s)
MA_Eve = 3.14; % mean anomaly (rad) at 0s UT

% Eve Moons 

% Gilly
mu_Gil = 8289449.8; % m^3/s^2
Tsr_Gil = 28255; % sidereal day (s)
Req_Gil = 13000; % equatorial radius (m)

% Gilly-Eve parameters
Ra_Gil = 48825000; % apoapsis (m)
Rp_Gil = 14175000; % periapsis (m)
e_Gil = 0.55; % orbital eccentricity (rad)
inc_Gil = 12*pi()/180; % orbital inclination (rad)
argp_Gil = 10*pi()/180; % argument of periapsis (rad)
lonasc_Gil = 80*pi()/180; % longitude of ascending node (rad)
Torb_Gil = 388587; % orbital period (s)
MA_Gil = 0.9; % mean anomaly (rad) at 0s UT

% Moho Constants  - - - - - - - - - - - - - - - - - - - - - 
mu_Moh = 1.6860938E11; % m^3/s^2
Tsr_Moh = 1210000; % sidereal day (s)
Req_Moh = 250000; % equatorial radius (m)

% solar parameters
Ra_Moh = 6315765981; % apoapsis (m)
Rp_Moh = 4210510628; % periapsis (m)
e_Moh = 0.2; % orbital eccentricity (rad)
inc_Moh = 7*pi()/180; % orbital inclination (rad)
argp_Moh = 15*pi()/180; % argument of periapsis (rad)
lonasc_Moh = 70*pi()/180; % longitude of ascending node (rad)
Torb_Moh = 2918346.4; % orbital period (s)
MA_Moh = 3.14; % mean anomaly (rad) at 0s UT

% Dres Constants - - - - - - - - - - - - - - - - - - - - - -
mu_Dre = 2.1484489E10; % m^3/s^2
Tsr_Dre = 34800; % sidereal day (s)
Req_Dre = 138000; % equatorial radius (m)

% solar parameters 
Ra_Dre = 46761053692; % apoapsis (m)
Rp_Dre = 34917642714; % periapsis (m)
e_Dre = 0.145;  % eccentricity 
inc_Dre = 5*pi()/180; % inclination (rad)
argp_Dre = 90*pi()/180; % argument of periapsis (rad)
lonasc_Dre = 280*pi()/180; % longitude of ascending node (rad)
Torb_Dre = 47893063;  % orbital period (s)
MA_Dre = 3.14;  % mean anomaly at 0s UT (rad)

% Jool Constants - - - - - - - - - - - - - - - - - - - - - -
mu_Joo = 2.8252800E14; % m^3/s^2
Tsr_Joo = 36000; % sidereal day (s)
Req_Joo = 6000000; % equatorial radius (m)

% solar parameters 
Ra_Joo = 72212238387; % apoapsis (m)
Rp_Joo = 65334882253; % periapsis (m)
e_Joo = 0.05; % eccentricity 
inc_Joo = 1.304*pi()/180; % inclination (deg)
argp_Joo = 0; % argument of periapsis (deg)
lonasc_Joo = 52*pi()/180; % longitude of ascending node (deg)
Torb_Joo = 104661432; % orbital period (s)
MA_Joo = 0.1; % mean anomaly at 0s UT 

% Jool Moons 

% Laythe
mu_Lay = 1.962E12; % m^3/s^2
Tsr_Lay = 52980.879; % sidereal day (s)
Req_Lay = 500000; % equatorial radius (m)

% Laythe-Jool Parameters
Ra_Lay = 27184000; % apoapsis (m)
Rp_Lay = 27184000; % periapsis (m)
e_Lay = 0;  % eccentricity 
inc_Lay = 0; % inclination (rad)
argp_Lay = 0; % argument of periapsis (rad)
lonasc_Lay = 0; % longitude of ascending node (rad)
Torb_Lay = 52981; % orbital period (s)
MA_Lay = 3.14; % mean anomaly at 0s UT 

% Vall
mu_Val = 2.0748150E11; % m^3/s^2
Tsr_Val = 105962.09; % s
Req_Val = 300000; 

% Vall-Jool Parameters 
Ra_Val = 43152000;
Rp_Val = 43152000;
e_Val = 0;
inc_Val = 0;
argp_Val = 0;
lonasc_Val = 0;
Torb_Val = 105962;
MA_Val = 0.9;

% Tylo
mu_Tyl = 2.8252800E12; 
Tsr_Tyl = 211926.36;
Req_Tyl = 600000;

% Tylo-Jool Parameters 
Ra_Tyl = 68500000;
Rp_Tyl = 68500000;
e_Tyl = 0;
inc_Tyl = 0.025;
lonasc_Tyl = 0;
Torb_Tyl = 212356.4;
MA_Tyl = 3.14;

% Bop
mu_Bop = 2.4868349E9;
Tsr_Bop = 544507.43;
Req_Bop = 65000;

% Bop-Jool Parameters 
Ra_Bop = 158697500;
Rp_Bop = 98302500;
e_Bop = 0.235;
inc_Bop = 15*pi()/180;
lonasc_Bop = 10*pi()/180;
Torb_Bop = 544507;
MA_Bop = 0.9;

% Pol
mu_Pol = 7.2170208E8;
Tsr_Pol = 901902.62;
Req_Pol = 44000;

% Pol-Jool Parameters 
Ra_Pol = 210624207;
Rp_Pol = 149155794;
e_Pol = 0.171
inc_Pol = 5.25*pi()/180;
lonasc_Pol = 2*pi()/180;
Torb_Pol = 901903;
MA_Pol = 0.9;















% Eeloo Constants - - - - - - - - - - - - - - - - - - - - - - 
mu_Elo = 7.4410815E10; % m^3/s^2
Tsr_Elo = 19460; % sidereal day (s)
Req_Elo = 210000; % equatorial radius (m)

% solar parameters 
Ra_Elo = 113549713200; % apoapsis (m)
Rp_Elo = 666879268900; % periapsis (m)
e_Elo = 0.26; % eccentricity 
inc_Elo = 6.15*pi()/180; % inclination (rad)
argp_Elo = 260*pi()/180; % argument of periapsis (rad)
lonasc_Elo = 50*pi()/180; % longitude of ascending node (rad)
Torb_Elo = 156992048; % orbital period (s)
MA_Elo = 3.14; % mean anomaly at 0s UT (rad) 



