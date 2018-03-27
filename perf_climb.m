%{
This function outputs the performance characteristics during the climb phase
---------------------------------------------------------------------------
Inputs:
alt_climb:         top of climb altitude [m]
CL_max:            max lift coefficient
W_climb:           Climb weight from [N]
Sref:              referrence area [m^2]
b_unswept:         unswept span [m]
---------------------------------------------------------------------------
Outputs:
ROC:               Rate of climb (m/s)
gamma_climb:       Climb angle (deg)
S_climb:           Distance covered during clim (m)
dt_climb:          Time to climb
===========================================================================
%}
function [ROC, gamma_climb, S_climb, dt_climb] = perf_climb(alt_climb, CL_max, W_climb, Sref, b_unswept)
%--------------------------------------------------------------------------
[~, ~, ~, rho_SL, son_SL, ~, ~, ~, ~, ~] = ATMO(0, 'm');           % sea level
[~, ~, ~, rho_top, son_top, ~, ~, ~, ~, ~] = ATMO(alt_climb, 'm'); % top of climb
rho_climb = 0.5*(rho_SL + rho_top);                                % everage air density [kg/m^3]
son_climb = 0.5*(son_SL + son_top);                                % everage speed of sound [m/s]
%--------------------------------------------------------------------------
CL_climb = CL_max/(1.25^2);	 % climb lift coefficient
gamma_climb = 15;            % [deg] climb angle

V_stall = sqrt(2*W_climb/(rho_climb*Sref*CL_climb)); % [m/s]
V_climb = sqrt((2/(rho_climb*CL_climb))*(W_climb/Sref)*cosd(gamma_climb)); % [Raymer eq. 17.38]

if V_climb <= V_stall
    V_climb = 1.5*V_stall; % [m/s]
end

M_climb = V_climb/son_climb; % climb Mach number
%--------------------------------------------------------------------------
% Calculate sweep for climb:
sweep_deg_climb = acosd(0.7/M_climb);
if sweep_deg_climb > 70   % Limit the sweep angle to 70 ish degrees
    sweep_deg_climb = 70;
end
b_swept_climb = b_unswept*cosd(sweep_deg_climb);  % Span at sweep angle [m]
AR_climb = b_swept_climb^2/Sref;                  % Swept aspect ratio
% ========================================================================
ROC = V_climb*sind(gamma_climb);                % rate of climb [m/s]
S_climb = alt_climb/tand(gamma_climb);          % [m] range covered to climb to subsonic cruise altitude
dt_climb = S_climb/(V_climb*cosd(gamma_climb)); % [s] time to climb to subsonic cruise altitude
%--------------------------------------------------------------------------
fprintf('\n\n ============================ Climb Results  ============================= \n');
fprintf('\n Climb sweep angle:              Gamma = %g [deg] ', sweep_deg_climb);
fprintf('\n Climb aspect ratio:             AR_climb = %g [deg] ', AR_climb);
fprintf('\n effective span during climb:    b_eff = %g [m] = %g [ft] ', b_swept_climb, b_swept_climb*3.2808399);
fprintf('\n\n -------------------------------------------------------------------- ');
fprintf('\n Rate of Climb:                  ROC = %g [m/s] = %g [ft/s]', ROC, ROC*3.2808399);
fprintf('\n Climb Angle:                    gamma = %g [deg] ', gamma_climb);
fprintf('\n Climb Velocity:                 V = %g [m/s] = %g [ft/s]', V_climb, V_climb*3.2808399);
fprintf('\n\n -------------------------------------------------------------------- ');
fprintf('\n Time to climb:                  TOF = %g [s] = %g [min]', dt_climb, dt_climb/60);
fprintf('\n\n ========================================================================= \n');
end
