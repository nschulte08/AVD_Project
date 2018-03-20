%{
This function returns the performance parameters during descent
---------------------------------------------------------------------------
Inputs:
alt_cr:        Altitude of cruise [m]
V_descend :    descent velocity [m/s]
---------------------------------------------------------------------------
Outputs:
S_descend:     Descent Distance Covered
dt_descend:    Time of descent phase
===========================================================================
%}
function [S_descend, dt_descend] = perf_descent(alt_cr, V_descend)

gamma_descend = -atand(1000/(3*5280));       % 3:1 rule
ROD = V_descend*sind(gamma_descend);         % [m/s] rate of descent
%--------------------------------------------------------------------------
S_descend = alt_cr/abs(tand(gamma_descend)); % [m] range covered during descent
dt_descend = S_descend/V_descend;            % [s] time to descend
%--------------------------------------------------------------------------
fprintf('\n\n ============================== Descent Results  ============================== \n');
fprintf('\n Descent from Subsonic Cruising Altitude: (h = %g [m])', alt_cr);
fprintf('\n Descent Angle:         gamma = %g [deg] ', gamma_descend);
fprintf('\n Rate of Descent:         ROD = %g [m/s] = %g [ft/s]', ROD, ROD*3.2808399);
fprintf('\n Descent Velocity:          V = %g [m/s] = %g [ft/s]', V_descend, V_descend*3.2808399);
fprintf('\n\n -------------------------------------------------------------------- ');
fprintf('\n Range Covered During Descent:   R = %g [km] = %g [miles]', S_descend/1000, S_descend*0.000621371);
fprintf('\n Time to Descend:               dt = %g [min]', dt_descend/60);
fprintf('\n\n -------------------------------------------------------------------- ');
end
