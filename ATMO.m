%{
===========================================================================
	STANDARD ATMOSPHERE

    by: Shawn McCullough
	last modified: 3/9/2017
%{
---------------------------------------------------------------------------
INPUTS:

h_g   = geometric altitude (dinstance above sea level)
        if units = 'M' , h_g must be in meters
        if units = 'E' , h_g must be in feet

units = string specifying which set of units are to be used
        acceptable inputs: 'M'  OR  'E' (not case sensitive)
        M = metric
        E = english (imperial)

display = string specifying if user wishes to display results to the
          command window. 'Y' = yes, (not case sensitive)
          any other input = no (a lack of input = no)

EXAMPLE:
[h, P, T, rho, son, mu, nu, DELTA, THETA, SIGMA] = ATMO(12500, 'M', 'Y') %display
[h, P, T, rho, son, mu, nu, DELTA, THETA, SIGMA] = ATMO(12500, 'E') % do not display
---------------------------------------------------------------------------
OUTPUTS:

h     = geopotential altitdue             [m]      OR  [ft]
P     = absolute pressure                 [Pa]     OR  [psi]
T     = absolute temperature              [K]      OR  [deg R]
rho   = density                           [kg/m^3] OR  [slug/ft^3]
son   = speed of sound                    [m/s]    OR  [ft/s]
mu    = absolute (dynamic) viscosity      [Pa*s]   OR  [psi*s]
nu    = kinematic viscosity               [m^2/s]  OR  [in^2/s]
DELTA = pressure ratio = P/P_reference
THETA = temperature ratio = T/T_reference
SIGMA = density ratio = rho/rho_reference
---------------------------------------------------------------------------
NOTES:
%{
This program assumes that the acceleration due to gravity is constant. This
assumption results in a simplified analysis, however it also results in an
approximate analysis, and the deviation from actual values will increase as
the altitude is increased. The specific heat ratio is also assumed to be a
constant, which further adds to the approximate nature of this program.
The location chosen for reference values is at 45 degrees latitude.
 
Equations and reference values were obtained from
"Introduction to Flight", 3rd Edition by John D. Anderson 
%}
%}
===========================================================================
%}
function [h, P, T, rho, son, mu, nu, DELTA, THETA, SIGMA] = ATMO(h_g, units, display)
%--------------------------------------------------------------------------
if nargin < 3 || isempty(display)
    display = 'N';
end
%--------------------------------------------------------------------------
if h_g < 0
    err_msg = sprintf(' Altitude must be greater than zero.');
	error('sytnax:requirements','\n\n %s \n\n',err_msg)
end
%--------------------------------------------------------------------------
if strcmp(units,'E')
    h_g = h_g*0.3048; % convert to meters for calculations
    % (all units will be converted back to english units at end of program)
end
%==========================================================================
    function [T, P, rho, son, mu, nu, DELTA, THETA, SIGMA] = gradient(a, T1, P1, rho1, h1, h)
        T = T1 + a*(h - h1);                                               % [K]
        P = P1*((T/T1)^(-g_0/(a*R)));                                      % [Pa]
        rho = rho1*((T/T1)^(-(g_0/(a*R) + 1)));                            % [kg/m^3]
        son = sqrt(gam*R*T);                                               % [m/s]
        mu = mu_0*((288.16 + 110)/(T + 110))*(T/288.16)^(3/2);             % [Pa*s] Sutherland's Law
        nu = mu/rho;                                                       % [m^2/s]
        DELTA = P/P_0;
        THETA = T/T_0;
        SIGMA = rho/rho_0;
    end
%==========================================================================
    function [T, P, rho, son, mu, nu, DELTA, THETA, SIGMA] = isothermal(T1, P1, rho1, h1, h)
        T = T1;
        P = P1*exp(-(g_0/(R*T))*(h - h1));                                 % [Pa]
        rho = rho1*exp(-(g_0/(R*T))*(h - h1));                             % [kg/m^3]
        son = sqrt(gam*R*T);                                               % [m/s]
        mu = mu_0*((288.16 + 110)/(T + 110))*(T/288.16)^(3/2);             % [Pa*s] Sutherland's Law
        nu = mu/rho;                                                       % [m^2/s]
        DELTA = P/P_0;
        THETA = T/T_0;
        SIGMA = rho/rho_0;
    end
%==========================================================================
% sea level conditions:
g_0 = 9.8;                 % [m/s^2]
T_0 = 288.16;              % [K]
P_0 = 101325;              % [Pa]
rho_0 = 1.225;             % [kg/m^3]
mu_0 = 1.7894e-5;          % [Pa*s]
R = 287;                   % [J/(kg*K)] gas constant
gam = 1.4;                 % specific heat ratio of air
%==========================================================================
R_m = 6.353766e6;          % [m] radius of the earth [m] (at 45 deg lat)
h = (R_m/(R_m + h_g))*h_g; % [m] geopotential altitude
%==========================================================================
    if h == 0
        T = T_0;
        P = P_0;
        rho = rho_0;
        son = sqrt(gam*R*T);
        mu = mu_0;
        nu = mu/rho;
        DELTA = 1.0;
        THETA = 1.0;
        SIGMA = 1.0;
    %----------------------------------------------------------------------
    elseif h > 0 && h < 11e3
        % initial conditions:
        T1 = T_0;          % [K]
        a = -6.5e-3;       % [K/m] lapse rate
        h1 = 0;            % [m]
        P1 = P_0;          % [Pa]
        rho1 = rho_0;      % [kg/m^3]
        %
        [T, P, rho, son, mu, nu, DELTA, THETA, SIGMA] = gradient(a, T1, P1, rho1, h1, h);
	%----------------------------------------------------------------------
    elseif h >= 11e3 && h < 25e3 %(isothermal layer)
        % initial conditions:
        T1 = 216.66;       % [K]
        h1 = 11e3;         % [m]
        P1 = 22650.2;      % [Pa]
        rho1 = 0.364205;   % [kg/m^3]
        %
        [T, P, rho, son, mu, nu, DELTA, THETA, SIGMA] = isothermal(T1, P1, rho1, h1, h);
    %----------------------------------------------------------------------
	elseif h >= 25e3 && h < 47e3
        % initial conditions:
        T1 = 216.66;       %[K]
        a = 3e-3;          %[K/m] lapse rate
        h1 = 25e3;         %[m]
        P1 = 2493.59;      %[Pa]
        rho1 = 0.0400957;  %[kg/m^3]
        %
        [T, P, rho, son, mu, nu, DELTA, THETA, SIGMA] = gradient(a, T1, P1, rho1, h1, h);
    %----------------------------------------------------------------------
	elseif h >= 47e3 && h < 53e3 %(isothermal layer)
        % initial conditions:
        T1 = 282.66;       % [K]
        h1 = 47e3;         % [m]
        P1 = 120.88;       % [Pa]
        rho1 = 0.00148985; % [kg/m^3]
        %
        [T, P, rho, son, mu, nu, DELTA, THETA, SIGMA] = isothermal(T1, P1, rho1, h1, h);
    %----------------------------------------------------------------------
	elseif h >= 53e3 && h < 79e3
        % initial conditions:
        T1 = 282.66;       %[K]
        a = -4.5e-3;       %[K/m] lapse rate
        h1 = 53e3;         %[m]
        P1 = 58.5556;      %[Pa]
        rho1 = 0.0007217;  %[kg/m^3]
        %
        [T, P, rho, son, mu, nu, DELTA, THETA, SIGMA] = gradient(a, T1, P1, rho1, h1, h);
    %----------------------------------------------------------------------
	elseif h >= 79e3 && h < 90e3 %(isothermal layer)
        % initial conditions:
        T1 = 165.66;        % [K]
        h1 = 79e3;          % [m]
        P1 = 1.01574;       % [Pa]
        rho1 = 2.13607e-05; % [kg/m^3]
        %
        [T, P, rho, son, mu, nu, DELTA, THETA, SIGMA] = isothermal(T1, P1, rho1, h1, h);
    %----------------------------------------------------------------------
	elseif h >= 90e3 && h <= 105e3
        % initial conditions:
        T1 = 165.66; %[K]
        a = 4e-3; %[K/m] lapse rate
        h1 = 90e3; %[m]
        P1 = 0.105216; %[Pa]
        rho1 = 2.21267e-06; %[kg/m^3]
        %
        [T, P, rho, son, mu, nu, DELTA, THETA, SIGMA] = gradient(a, T1, P1, rho1, h1, h);
    %----------------------------------------------------------------------
    elseif h > 105e3
        err_msg = sprintf(' Altitude exceeds maximum value. \n For this program, the altitude must be less than 105 km (344488.19 ft)');
        error('sytnax:requirements','\n\n %s \n\n',err_msg)
    %----------------------------------------------------------------------
    end
%==========================================================================
if strcmpi(units,'E')
    % convert to english units:
    h_g = h_g*3.2808399;    % [ft]
    h = h*3.2808399;        % [ft];
    P = P*0.00014503762;    % [psi]
    T = T*1.8;              % [deg R]
    rho = rho*0.0019403196; % [slug/ft^3]
    son = son*3.2808399;    % [ft/s]
    mu = mu*0.00014503762;  % [psi*s]
    nu = nu*1549.9959;      % [in^2/s]
    
    if strcmpi(display,'Y')
        fprintf('\n\n =============================================== \n');
        fprintf('\n At a geometric altitude of h_g = %g [ft] \n', h_g);
        fprintf('\n Geopotential altitude:  h = %g [ft] \n', h);
        fprintf('\n Temperature:         T   = %g [deg R] ', T);
        fprintf('\n Pressure:            P   = %g [psi] ', P);
        fprintf('\n Density:             rho = %g [slug/ft^3] ', rho);
        fprintf('\n Speed of Sound:      a   = %g [ft/s] ', son);
        fprintf('\n Absolute Viscosity:  mu  = %g [psi*s] ', mu);
        fprintf('\n Kinematic Viscosity: nu  = %g [in^2/s] \n', nu);
        fprintf('\n Temperature Ratio:   THETA = %g', THETA);
        fprintf('\n Pressure Ratio:      DELTA = %g', DELTA);
        fprintf('\n Density Ratio:       SIGMA = %g', SIGMA);
        fprintf('\n\n =============================================== \n\n');
    end
%--------------------------------------------------------------------------    
else %(metric units):
    if strcmpi(display,'Y')
        fprintf('\n\n =============================================== \n');
        fprintf('\n At a geometric altitude of h_g = %g [km] \n', h_g/1000);
        fprintf('\n Geopotential altitude:  h = %g [km] \n', h/1000);
        fprintf('\n Temperature:         T   = %g [K] ', T);
        fprintf('\n Pressure:            P   = %g [kPa] ', P/1000);
        fprintf('\n Density:             rho = %g [kg/m^3] ', rho);
        fprintf('\n Speed of Sound:      a   = %g [m/s] ', son);
        fprintf('\n Absolute Viscosity:  mu  = %g [Pa*s] ', mu);
        fprintf('\n Kinematic Viscosity: nu  = %g [m^2/s] \n', nu);
        fprintf('\n Temperature Ratio:   THETA = %g', THETA);
        fprintf('\n Pressure Ratio:      DELTA = %g', DELTA);
        fprintf('\n Density Ratio:       SIGMA = %g', SIGMA);
        fprintf('\n\n =============================================== \n\n');
    end
end
%==========================================================================
end