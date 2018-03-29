%{
This function returns the performance parameters during descent
---------------------------------------------------------------------------
Inputs:
alt_cr:        Altitude of cruise [m]
V_descend :    descent velocity [m/s]
M_perp:        wing perpendicular Mach number
---------------------------------------------------------------------------
Outputs:
S_descend:     Descent Distance Covered
dt_descend:    Time of descent phase
sweep_descend: Sweep schedule = [sweep_deg_top, sweep_deg_bottom];
===========================================================================
%}
function [S_descend, dt_descend, sweep_descend] = perf_descent(alt_cr, M_perp)
%--------------------------------------------------------------------------
[~, ~, ~, ~, son_SL, ~, ~, ~, ~, ~] = ATMO(0, 'm');       % sea level
[~, ~, ~, ~, son_top, ~, ~, ~, ~, ~] = ATMO(alt_cr, 'm'); % top of descent

M = 0.78; % maximum Mach number
V_descend = M*son_top; % [m/s]
%--------------------------------------------------------------------------
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
%% ========================================================================
% Calculate sweep schedule for descent:
%--------------------------------------------------------------------------
% top of descent:
M_top = V_descend/son_top; % top of climb Mach number
if M_perp < M_top % to avoid complex numbers
    sweep_deg_top = acosd(M_perp/M_top); % This has to be limited to 70 ish degrees!
    if sweep_deg_top > 70                % Limit the sweep angle
        sweep_deg_top = 70;
    end
else
    sweep_deg_top = 0;
end
%--------------------------------------------------------------------------
% bottom of descent:
M_bottom = V_descend/son_SL; % bottom of climb Mach number
if M_perp < M_bottom % to avoid complex numbers
    sweep_deg_bottom = acosd(M_perp/M_bottom); % This has to be limited to 70 ish degrees!
    if sweep_deg_bottom > 70                % Limit the sweep angle
        sweep_deg_bottom = 70;
    end
else
    sweep_deg_bottom = 0;
end
%--------------------------------------------------------------------------
sweep_descend = [sweep_deg_top, sweep_deg_bottom];
%--------------------------------------------------------------------------
end
