% Damon Printz
% 4/2/2020
% calculates time until Kerbal Space Center is at the current 
% value of Right Ascension

% Input target right ascension RA in degrees 
% Input current time t (curTime = [year, day, hour, minute, second])

% get target time T ( [year, day, hour, minute, second] )
% get time until target t_until in seconds
% current right ascension

function [T, t_until, RAcur] = get_KSC_t_until_ra (RA, curTime)

% kerbal space center longitude in degrees, West = (-)
ksp_lon = -(74+34/60+31/3600);

Tsidereal = 21549.425; % sidereal day (s)
Tday = 21600; % solar day (s)
Tyear = 9203545; % Kerbin year (s)

% calculate current time 
tc = curTime(1)*Tyear+curTime(2)*Tday+curTime(3)*3600+curTime(4)*60+curTime(5);

% get current RA 
RAcur = mod(tc, Tsidereal) * 360 / Tsidereal + ksp_lon;

% calculate the change in time required to reach the destination right ascension 
% in seconds 
t_until = ( RA - RAcur ) * ( Tsidereal / 360 );

% if the calculator says to go backwards in time, simply add one orbital period 
% to the mix 
while(t_until < 0)
t_until = t_until + Tsidereal;
endwhile 

% calculate target time and convert to UT 
tt = tc + t_until;

T = zeros(size(curTime));

T(1) = floor( tt / Tyear ); % get years 
T(2) = floor( ( tt - T(1) * Tyear ) / Tday); % get days 
T(3) = floor( ( tt - T(1) * Tyear - T(2) * Tday ) / 3600 ); % get hours  
T(4) = floor( ( tt - T(1) * Tyear - T(2) * Tday - T(3) * 3600 ) / 60 ); % get minutes 
T(5) = tt - T(1) * Tyear - T(2) * Tday - T(3) * 3600 - T(4) * 60; % get seconds 

% print target universal time 
fprintf('Target UT: y%.4i:d%.3i - h%i:m%i:s%.2i\n',T(1),T(2),T(3),T(4),T(5))
% print time until target intersect 
fprintf('T until target: %.1f s\n', t_until)
% print current right ascension 
fprintf('Current RA: %.1f deg\n', RAcur)



endfunction
