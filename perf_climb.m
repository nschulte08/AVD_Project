%{
This function outputs the performance characteristics during the climb phase
---------------------------------------------------------------------------
Inputs:
CL_climb:          Climb lift coefficient
CD_climb:          Climb drag coefficient
V_stall:           Stall velocity (m/s)
V_sub_cr:           Subsonic cruise velocity (m/s)
rho_TO:            Takeoff air density (kg/m^3)
T_max:             Sea level max thrust (N)
W_climb:           Climb weight from wt fraction analysis
---------------------------------------------------------------------------
Outputs:
ROC:               Rate of climb (m/s)
gamma_climb:       Climb angle (deg)
S_climb:           Distance covered during clim (m)
dt_climb:          Time to climb
===========================================================================
%}
function [ROC, gamma_climb, S_climb, dt_climb] = perf_climb(CL_climb, CD_climb, V_stall, V_sub_cr, rho_TO, T_max, W_climb, alt_climb, Sref)


%% ========================================================================
% Maximum ROC:
%--------------------------------------------------------------------------
V = linspace(V_stall,V_sub_cr)';                   % [m/s] veloctiy vector
D = 0.5*rho_TO*V.^2*CD_climb*Sref;                 % [N] drag
ROC = V.*(T_max - D)./W_climb;                     % [m/s] rate of climb
[max_ROC, ind] = max(ROC);                         % [m/s] maximum rate of climb

gamma_climb_max = asind((T_max - D(ind))/W_climb); % [deg] corresponding climb angle
gamma_climb = floor(gamma_climb_max);              % [deg] climb angle
%gamma_climb = 15; % [deg]

%{
figure('Name','Max ROC vs Velocity at Sea Level','NumberTitle','off');
plot(V,ROC,'k', 'LineWidth',2);
title('Maximum Rate of Climb','FontSize',18);
xlabel('Airspeed [m/s]','FontSize',12);
ylabel('ROC [m/s]','FontSize',12);
%}
% =========================================================================
fprintf('\n\n ============================ Climb Results  ============================= \n');
fprintf('\n Sea level Rate of Climb: \n');
fprintf('\n Maximum Rate of Climb:         ROC = %g [m/s] = %g [ft/s]', max_ROC, max_ROC*3.2808399);
fprintf('\n Corresponding Climb Angle:     gamma = %g [deg] ', gamma_climb_max);
fprintf('\n Corresponding Climb Velocity:  V = %g [m/s] = %g [ft/s]', V(ind), V(ind)*3.2808399);
fprintf('\n\n -------------------------------------------------------------------- ');
%--------------------------------------------------------------------------
% Normal Climb:
V_climb = sqrt((2/(rho_TO*CL_climb))*(W_climb/Sref)*cosd(gamma_climb)); % [Raymer eq. 17.38]

% ROC_norm = V_climb*sind(gamma_climb); % [m/s]

S_climb = alt_climb/tand(gamma_climb);          % [m] range covered to climb to subsonic cruise altitude
dt_climb = S_climb/(V_climb*cosd(gamma_climb)); % [s] time to climb to subsonic cruise altitude

end

