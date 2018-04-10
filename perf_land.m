%{
This function outputs the performance during landing
---------------------------------------------------------------------------
Inputs:
alt_land:   altitude of landing runway [m]
Sref:       referrence area [m^2]
AR_low:     aspect ratio
W_land:     landing weight [N]  
CL_max:     max lift coefficient
T_max:      max thrust
TR:         taper ratio
M_land:     landing Mach number
---------------------------------------------------------------------------
Outputs:
S_land:     landing distance [m]
FAR_land:   FAR landing distance [m]
===========================================================================
%}
function [S_land, FAR_land, V_TD] = perf_land(alt_land, Sref, AR_low, W_Land, CL_max, TR, SM)

[~, ~, ~, rho_land, son_land, ~, ~, ~, ~, ~] = ATMO(alt_land, 'M');

mu = 0.5;                 % friction coefficient (Yechout p.108)
g = 9.81;                 % [m/s^2] gravity
h_obst = 35*0.3048;       % [m] 35ft obstacle
gamma_land = 3;           % [deg] approach angle

V_stall = sqrt(2*W_Land/(rho_land*Sref*CL_max)); % [m/s] stall speed
V_approach = 1.2*V_stall;                        % [m/s] approach velocity
V_TD = 1.15*V_stall;                             % [m/s] touch down velocity

M_land = V_approach/son_land; % landing Mach number

CL_Land = W_Land/(Sref*0.5*rho_land*V_approach^2);
[CD_Land, CD_0_land] = aerofunk_drag_2(alt_land, M_land, Sref, CL_Land, SM, AR_low, TR);
%--------------------------------------------------------------------------
% flare:
V_flare = (V_approach + V_TD)/2;                  % [m/s]
n = 1.2;                                          % load factor
Radius_flare = V_flare^2/(g*(n-1));               % [m] flare radius
h_flare = Radius_flare*(1 - cosd(gamma_land));    % [m] flare height

S_flare = Radius_flare*sind(gamma_land);          % [m] distance covered during approach
%--------------------------------------------------------------------------
% approach:
S_approach = (h_obst - h_flare)/tand(gamma_land); % [m] distance covered during flare
%--------------------------------------------------------------------------
% free roll:
S_Freeroll = 3*V_TD; % [m] free roll distance
%--------------------------------------------------------------------------
% breaking distance:
S_brake = (W_Land/(g*rho_land*Sref*(CD_Land-mu*CL_Land)))*log(1 + (rho_land*Sref/(2*W_Land))*(CD_Land/mu - CL_Land)*V_TD^2); % Nicolai
%{
S_brake = (W_land/(g*rho_land*Sref*(CD_Land-mu*CL_Land)))*log(1 + ((CD_Land-mu*CL_Land)*rho_land*Sref*V_TD^2)/(2*(-T_Land + mu*W_land))); % Yechout
if imag(S_brake) ~= 0 || real(S_brake) < 0
    S_brake = 0;
end
%}
%--------------------------------------------------------------------------
% total landing length:
S_land = S_approach + S_flare + S_Freeroll + S_brake; % [m]
FAR_land = S_land*(1 + 2/3); % FAR requirement
%--------------------------------------------------------------------------
fprintf('\n\n ============================== Landing Results  ============================== \n');
fprintf('\n Approach Angle:      gamma = %g [deg]', gamma_land);
fprintf('\n Approach Speed:      V_a   = %g [m/s] = %g [ft/s]', V_approach, V_approach*3.2808399);
fprintf('\n Flare Speed:         V_f   = %g [m/s] = %g [ft/s]', V_flare, V_flare*3.2808399);
fprintf('\n Touch-Down Speed:    V_TD  = %g [m/s] = %g [ft/s]', V_TD, V_TD*3.2808399);
fprintf('\n\n ----------------------------------------------------------------------------- \n');
fprintf('\n Flare height:    h_flare = %g [m] = %g [ft]', h_flare, h_flare*3.2808399);
fprintf('\n Obstacle height: h_obst  = %g [m] = %g [ft]', h_obst, h_obst*3.2808399);
fprintf('\n Flare radius:    R = %g [m] = %g [ft]', Radius_flare, Radius_flare*3.2808399);
fprintf('\n\n ----------------------------------------------------------------------------- \n');
fprintf('\n Distance covered on approach:    S_Approach  = %g [m] = %g [ft]', S_approach, S_approach*3.2808399);
fprintf('\n Distance covered during flare:   S_Flare     = %g [m] = %g [ft]', S_flare, S_flare*3.2808399);
fprintf('\n Free Roll Distance:              S_Free_roll = %g [m] = %g [ft]', S_Freeroll, S_Freeroll*3.2808399);
fprintf('\n Breaking Distance:               S_Breaking  = %g [m] = %g [ft]', S_brake, S_brake*3.2808399);
fprintf('\n\n Total Landing Distance:        S_Land = %g [m] = %g [ft]', S_land, S_land*3.2808399);
fprintf('\n\n FAR Landing Distance:          S_FAR  = %g [m] = %g [ft]', FAR_land, FAR_land*3.2808399);
fprintf('\n\n =============================================================================== \n');
end
