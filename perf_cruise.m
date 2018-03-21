%{
This function outputs the performance during cruise
---------------------------------------------------------------------------
Inputs:
M_cr:         Cruise Mach
alt_cr:       Altitude for cruise
W_cr_start:   Weight at beginning of cruise
W_cr_end:     Weight at end of cruise
Sref:         Reference wing area (m^2)
SM:           static margin
AR:           Aspect Ratio
TSFC:         TSFC from propulsion
---------------------------------------------------------------------------
Outputs:
R_constH:     Constant altitude range (m)
R_CC:         Cruise climb range (m)
TOF_constH    Constant altitude time of flight
TOF_CC        Cruise climb time of flight
===========================================================================
%}
function [R_constH, R_CC, TOF_constH, TOF_CC] = perf_cruise(M_cr, alt_cr, W_cr_start, W_cr_end, Sref, SM, AR, TSFC, TR)

[~, ~, ~, rho_cr, son_cr, ~, ~, ~, ~, ~] = ATMO(alt_cr, 'M');
V_cr = M_cr*son_cr; % [m/s]

CL_cr = W_cr_start/(Sref*0.5*rho_cr*V_cr^2);     % lift coefficient cruise
[CD_cr, ~] = aerofunk_drag_2(alt_cr, M_cr, Sref, CL_cr, SM, AR, TR);

%% ========================================================================
% Constant altitude cruise:
R_constH = sqrt(2/(rho_cr*Sref))*(2/TSFC)*(sqrt(CL_cr)/CD_cr)*(sqrt(W_cr_start) - sqrt(W_cr_end)); % [m] Breguet range equation
%--------------------------------------------------------------------------
% Cruise climb:
R_CC = sqrt(2*W_cr_start/(rho_cr*Sref))*(1/TSFC)*(sqrt(CL_cr)/CD_cr)*log(W_cr_start/W_cr_end); % [m]
%--------------------------------------------------------------------------
% Time of flight:
TOF_constH = R_constH/V_cr; % [s]
TOF_CC = R_CC/V_cr;         % [s]
%% ========================================================================
fprintf('\n\n ============================== Cruise Results  =============================== \n');
if M_cr < 1
    fprintf('\n Subsonic cruise: \n');
else
    fprintf('\n Supersonic cruise: \n');
end
    fprintf('\n Cruise Mach number: M  = %g', M_cr); 
    fprintf('\n Cruise Velocity:    V  = %g [m/s] = %g [ft/s]', V_cr, V_cr*3.2808399);
    fprintf('\n Cruise Altitude:    h  = %g [m] = %g [ft] \n', alt_cr, alt_cr*3.2808399);
    fprintf('\n Constant Altitude Range:    R  = %g [km]   = %g [miles]', R_constH/1000, R_constH*0.000621371);
    fprintf('\n Cruise Climb Range:         R  = %g [km]   = %g [miles] \n', R_CC/1000, R_CC*0.000621371);
    fprintf('\n Constant Altitude Cruise Time of Flight:   TOF  = %g [min]', TOF_constH/60);
    fprintf('\n Cruise Climb Time of Flight:               TOF  = %g [min]', TOF_CC/60);
    fprintf('\n\n ============================================================================== \n');
