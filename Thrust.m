%{
===========================================================================
	MAE 4350 - Aerospace Senior Design Capstone Project
    University of Texas at Arlington
    AS2 Supersonic Business Jet
---------------------------------------------------------------------------
    Thrust at a given Altitude and Mach Number
    ------------------------------------------
    by: Shawn McCullough
	last modified: 10/18/2017
%{
---------------------------------------------------------------------------
INPUTS:
h = altitude [m]
M = Mach number
T_max = max thrust [N]
---------------------------------------------------------------------------
OUTPUTS:
T = thrust [N]
---------------------------------------------------------------------------
%}
===========================================================================
%}
function [T] = Thrust(h, M, T_max)
%--------------------------------------------------------------------------
%T_max = 12345;     % [N] max installed thrust (total with all engines)
Theta_break = 1;   % theta break
TR = Theta_break;  % throttle ratio
%--------------------------------------------------------------------------
[~, ~, ~, ~, ~, ~, ~, DELTA, THETA, ~] = ATMO(h, 'M'); % atmospheric properties
THETA_0 = THETA*(1+((1.4-1)/2)*M^2);                   % temperature ratio
DELTA_0 = DELTA*(1+((1.4-1)/2)*M^2)^(1.4/(1.4-1));     % pressure ratio
%--------------------------------------------------------------------------
 % Mattingly eq. 2.54.a :
if THETA_0 <= TR
    Alpha = DELTA_0;
else
    Alpha = DELTA_0*(1 - 3.5*(THETA_0-TR)/THETA_0);
end
%--------------------------------------------------------------------------
T = Alpha*T_max; % [N] thrust
%--------------------------------------------------------------------------
end
