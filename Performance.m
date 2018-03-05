%{
===========================================================================
	MAE 4351 - Aerospace Senior Design Capstone Project
    University of Texas at Arlington
---------------------------------------------------------------------------
    Performance Calculations
    ------------------------
    by: Shawn McCullough
	last modified: 02/28/2017
===========================================================================
INPUTS:
% Known flight parameters:
h_sub_R = 13106;   % [m] subsonic cruise altitude
M_sub_R = 0.95;    % subsonic cruise Mach number
h_super_R = 15545; % [m] supersonic cruise altitude
M_super_R = 1.4;   % supersonic cruise Mach number
%--------------------------------------------------------------------------
% geometry:
S_w = 125;         % [m^2] refernce area (wing planform)
b = 18.6;          % [m] wing span
%--------------------------------------------------------------------------
% Propulsion
n = 4; % number of engines
TSFC_sub   = 0.###/3600; % [1/s]
TSFC_super = 0.###/3600; % [1/s]
T_max = #######;         % [N] max installed thrust
%--------------------------------------------------------------------------
% Aerodynamics:
CL_max = 1.2;        % max lift coefficient
CL_Land = 0.745;     % landing lift coefficient
e = 0.8;             % Oswald efficiency factor (assumption)
%--------------------------------------------------------------------------
% weights:
weights_sub   = structure of weights for subsonic mission
weights_super = structure of weights for supersonic mission
===========================================================================
%}
function Performance(h_sub_R, M_sub_R, h_super_R, M_super_R, S_w, b, e, CL_max, CL_Land, T_max, TSFC_super, TSFC_sub, n, weights_sub, weights_super)
%% ========================================================================
g = 9.81;        % [m/s^2] gravity
AR = b^2/S_w;    % apsect ratio
K = 1/(pi*e*AR); % drag thingy
%--------------------------------------------------------------------------
% weights (from Ben's weight script turned into simplified function):
MTOW = 4.44822*weights_sub.W_to;                                               % [N] max take off weight
W_climb_sub = 4.44822*(weights_sub.W_TO_end + weights_sub.W_climb_end)/2;      % [N] average value
%W_climb_sup = 4.44822*(weights_super.W_TO_end + weights_super.W_climb_end)/2; % [N] average value
W_climb = W_climb_sub;                                                         % [N] average value
W_cruise_end = 4.44822*weights_sub.W_cruise_end;                               % [N] weight at end of subsonic cruise
W_cruise_start_sub = 4.44822*weights_sub.W_climb_end;                          % [N] weight at start of subsonic cruise
W_cruise_start_sup = 4.44822*weights_super.W_accel_end;                        % [N] weight at start of supersonic cruise (end of acceleration phase)
%W_descend = 4.44822*(weights_sub.W_cruise_end + weights_sub.W_descent_end)/2; % [N] average value
W_land = 4.44822*weights_sub.W_land_end;                                       % [N] landing weight

%% ========================================================================
% Take off:
%--------------------------------------------------------------------------
T_TO = ((n-1)/n)*T_max; % [N] OEI scenario
% ground roll:
mu = (0.02 + 0.3)/2;                                 % friction coefficient [average value] (Yechout p.99)
[~, ~, ~, rho_TO, son_TO, ~, ~, ~, ~, SIGMA_TO] = ATMO(0, 'M');
V_stall = sqrt(2*MTOW/(rho_TO*S_w*CL_max));          % [m/s]
V_TO = 86;                                           % [m/s]
M_TO = V_TO/son_TO;                                  % Mach number
CD_0 = CD_0_funk(M_TO);                              % parasitic drag coefficient
q_TO = 0.5*rho_TO*V_TO^2;                            % [N/m^2] dynamic pressure
CL_opt = mu/(2*K);                                   % optimum lift coefficient for T.O.
CD_avg = CD_0 + K*CL_opt^2;                          % example 3.5 from Yechout (p.103)
D_avg = CD_avg*S_w*q_TO;                             % [N] average drag
L_avg = CL_opt*S_w*q_TO;                             % [N] average lift
F_r = mu*(MTOW - L_avg);                             % [N] rolling resistance
a_TO = (g/MTOW)*(T_TO - D_avg - F_r);                % [m/s^2] average acceleration
S_G = V_TO^2/(2*a_TO);                               % [m] ground distance
%--------------------------------------------------------------------------
% rotation:
t_R = 2;        % [s] rotation time
S_R = t_R*V_TO; % [m] distance covered during rotation
%--------------------------------------------------------------------------
% transition:
n = 1.2;                               % load factor
V_TR = 1.5*V_stall;                    % [m/s] transition velocity
Radius_TO = V_TR^2/(g*(n-1));          % [m] radius of transition
gamma_c = asind((T_TO - D_avg)/MTOW);  % [deg] climb angle
S_TR = Radius_TO*sind(gamma_c);        % [m] transition distance
h_TR = Radius_TO*(1 - cosd(gamma_c));  % [m] altitude gained during transition
%--------------------------------------------------------------------------
% climb:
h_obst = 35*0.3048; % [m] 35ft obstacle
if h_TR > h_obst
    S_C = 0; % [m] distance covered during climb
else
    S_C = (h_obst - h_TR)/tand(gamma_c);   % [m] distance covered during climb
end
%--------------------------------------------------------------------------
% total take-off distance:
S_TO = S_G + S_R + S_TR + S_C; % [m]
%--------------------------------------------------------------------------
% balanced field length: (Raymer eq. 17.110)
gamma_min = 0.027; % [deg]
G = (gamma_c - gamma_min)*pi/180;
U = 0.01*CL_max + 0.02;
BFL = (0.863/(1+2.3*G))*((MTOW/S_w)/(rho_TO*g*CL_opt) + h_obst)*(1/((T_TO/MTOW) - U) + 2.7*0.3048) + (655*0.3048/sqrt(SIGMA_TO)); % [m] balanced field length
%--------------------------------------------------------------------------
fprintf('\n\n =========================== Take-off Results  =========================== \n');
fprintf('\n Max Take-off Gross Weight: MTOW = %g [N] = %g [lbf] \n', MTOW, MTOW*0.22480902);
fprintf('\n Stall Speed:    V_s  = %g [m/s] = %g [ft/s]', V_stall, V_stall*3.2808399);
fprintf('\n Take-off Speed: V_TO = %g [m/s] = %g [ft/s]', V_TO, V_TO*3.2808399);
fprintf('\n Obstacle clearance height: h_obst = %g [m] = %g [ft]', h_obst, h_obst*3.2808399);
fprintf('\n\n ------------------------------------------------------------- \n');
fprintf('\n Ground Roll:         S_G  = %g [m] = %g [ft]', S_G, S_G*3.2808399);
fprintf('\n Rotation Distance:   S_R  = %g [m] = %g [ft]', S_R, S_R*3.2808399);
fprintf('\n Transition Distance: S_TR = %g [m] = %g [ft]', S_TR, S_TR*3.2808399);
fprintf('\n Climb Distance:      S_C  = %g [m] = %g [ft]', S_C, S_C*3.2808399);
fprintf('\n\n Total Take-off Distance:      S_TO  = %g [m] = %g [ft]', S_TO, S_TO*3.2808399);
fprintf('\n\n Balanced Field Length:         BFL  = %g [m] = %g [ft]', BFL, BFL*3.2808399);
fprintf('\n\n ========================================================================= \n');
%% ========================================================================
% Subsonic Cruise:
%--------------------------------------------------------------------------
%TSFC_sub = TSFC_sub*g; % [1/s] required form for range equations
%--------------------------------------------------------------------------
% Constant altitude cruise:
[~, ~, ~, rho_sub_R, son_sub_R, ~, ~, ~, ~, ~] = ATMO(h_sub_R, 'M');
V_sub_R = M_sub_R*son_sub_R; % [m/s]

CL_R_sub = W_cruise_start_sub/(S_w*0.5*rho_sub_R*V_sub_R^2);     % lift coefficient for subsonic cruise (R=range)
CD_R_sub = Drag_Buildup_FUNKY(CL_R_sub, 0.95);
%CD_R_sub = 0.01571;    % drag coefficient for subsonic cruise
%CL_R_sub = 0.0586;     % lift coefficient for subsonic cruise (R=range)

R_constH_sub = sqrt(2/(rho_sub_R*S_w))*(2/TSFC_sub)*(sqrt(CL_R_sub)/CD_R_sub)*(sqrt(W_cruise_start_sub) - sqrt(W_cruise_end)); % [m] Breguet range equation
%--------------------------------------------------------------------------
% Cruise climb:
R_CC_sub = sqrt(2*W_cruise_start_sub/(rho_sub_R*S_w))*(1/TSFC_sub)*(sqrt(CL_R_sub)/CD_R_sub)*log(W_cruise_start_sub/W_cruise_end); % [m]
%--------------------------------------------------------------------------
% Time of flight:
TOF_constH_sub = R_constH_sub/V_sub_R; % [s]
TOF_CC_sub = R_CC_sub/V_sub_R;         % [s]

%% ========================================================================
% Supersonic Cruise:
%--------------------------------------------------------------------------
%TSFC_super = TSFC_super*g; % [1/s] required form for range equations
%--------------------------------------------------------------------------
% Constant altitude cruise:
[~, ~, ~, rho_super_R, son_super_R, ~, ~, ~, ~, ~] = ATMO(h_super_R, 'M');
V_super_R = M_super_R*son_super_R; % [m/s]

CL_R_super = W_cruise_start_sup/(S_w*0.5*rho_super_R*V_super_R^2);  % lift coefficient for supersonic cruise
CD_R_super = Drag_Buildup_FUNKY(CL_R_super, 1.4);
%CL_R_super = 0.02698;  % lift coefficient for supersonic cruise
%CD_R_super = 0.01723;  % drag coefficient for supersonic cruise

R_constH_super = sqrt(2/(rho_super_R*S_w))*(2/TSFC_super)*(sqrt(CL_R_super)/CD_R_super)*(sqrt(W_cruise_start_sup) - sqrt(W_cruise_end)); % [m] Breguet range equation
%--------------------------------------------------------------------------
% Cruise climb:
R_CC_super = sqrt(2*W_cruise_start_sup/(rho_super_R*S_w))*(1/TSFC_super)*(sqrt(CL_R_super)/CD_R_super)*log(W_cruise_start_sup/W_cruise_end); % [m]
%--------------------------------------------------------------------------
% Time of flight:
TOF_constH_super = R_constH_super/V_super_R; % [s]
TOF_CC_super = R_CC_super/V_super_R;         % [s]

%% ========================================================================
% Maximum ROC:
%--------------------------------------------------------------------------
CL_climb = CL_max/(1.25^2); % climb lift coefficient
CD_climb = Drag_Buildup_FUNKY(CL_climb, 0.8);

V = linspace(V_stall,V_sub_R)';                    % [m/s] veloctiy vector
D = 0.5*rho_TO*V.^2*CD_climb*S_w;                  % [N] drag
ROC = V.*(T_max - D)./W_climb;                     % [m/s] rate of climb
[max_ROC, ind] = max(ROC);                         % [m/s] maximum rate of climb

V_climb = floor(V(ind));                           % [m/s] climb velocity
M_climb = V_climb/son_TO;                          % climb Mach number
gamma_climb_max = asind((T_max - D(ind))/W_climb); % [deg] corresponding climb angle
gamma_climb = floor(gamma_climb_max);              % [deg] climb angle
%gamma_climb = 15; % [deg]

figure('Name','Max ROC vs Velocity at Sea Level','NumberTitle','off');
plot(V,ROC,'k', 'LineWidth',2);
title('Maximum Rate of Climb','FontSize',18);
xlabel('Airspeed [m/s]','FontSize',12);
ylabel('ROC [m/s]','FontSize',12);
% =========================================================================
fprintf('\n\n ============================ Climb Results  ============================= \n');
fprintf('\n Sea level Rate of Climb: \n');
fprintf('\n Maximum Rate of Climb:         ROC = %g [m/s] = %g [ft/s]', max_ROC, max_ROC*3.2808399);
fprintf('\n Corresponding Climb Angle:     gamma = %g [deg] ', gamma_climb_max);
fprintf('\n Corresponding Climb Velocity:  V = %g [m/s] = %g [ft/s]', V(ind), V(ind)*3.2808399);
fprintf('\n\n -------------------------------------------------------------------- ');
%--------------------------------------------------------------------------
% Normal Climb:
%V_climb = sqrt((2/(rho_TO*CL_climb))*(W_climb/S_w)*cosd(gamma_climb)); % [Raymer eq. 17.38]

ROC_norm = V_climb*sind(gamma_climb); % [m/s]

S_climb_sub = h_sub_R/tand(gamma_climb);     % [m] range covered to climb to subsonic cruise altitude
S_climb_super = h_super_R/tand(gamma_climb); % [m] range covered to climb to supersonic cruise altitude

dt_climb_sub = S_climb_sub/(V_climb*cosd(gamma_climb));     % [s] time to climb to subsonic cruise altitude
dt_climb_super = S_climb_super/(V_climb*cosd(gamma_climb)); % [s] time to climb to supersonic cruise altitude
%--------------------------------------------------------------------------
fprintf('\n -------------------------------------------------------------------- ');
fprintf('\n For a Climb Angle of: gamma = %g [deg] \n', gamma_climb);
fprintf('\n ROC = %g [m/s] = %g [ft/s]', ROC_norm, ROC_norm*3.2808399);
fprintf('\n V   = %g [m/s] = %g [ft/s]', V_climb, V_climb*3.2808399);
fprintf('\n -------------------------------------------------------------------- ');
fprintf('\n Time to Reach Subsonic Cruising Altitude: (%g [m])', h_sub_R);
fprintf('\n dt = %g [min]', dt_climb_sub/60);
fprintf('\n Range Covered During Climb:\n R = %g [km] = %g [miles]', S_climb_sub/1000, S_climb_sub*0.000621371);
fprintf('\n -------------------------------------------------------------------- ');
fprintf('\n Time to Reach Supersonic Cruising Altitude: (%g [m])', h_super_R);
fprintf('\n dt = %g [min]', dt_climb_super/60);
fprintf('\n Range Covered During Climb:\n R = %g [km] = %g [miles]', S_climb_super/1000, S_climb_super*0.000621371);
fprintf('\n\n ============================================================================== \n');
%% ========================================================================
% cruise results:
fprintf('\n\n ============================== Cruise Results  =============================== \n');
fprintf('\n Subsonic cruise: \n');
fprintf('\n Cruise Mach number: M  = %g', M_sub_R); 
fprintf('\n Cruise Velocity:    V  = %g [m/s] = %g [ft/s]', V_sub_R, V_sub_R*3.2808399);
fprintf('\n Cruise Altitude:    h  = %g [m] = %g [ft] \n', h_sub_R, h_sub_R*3.2808399);
fprintf('\n Constant Altitude Range:    R  = %g [km]   = %g [miles]', R_constH_sub/1000, R_constH_sub*0.000621371);
fprintf('\n Cruise Climb Range:         R  = %g [km]   = %g [miles] \n', R_CC_sub/1000, R_CC_sub*0.000621371);
fprintf('\n Constant Altitude Cruise Time of Flight:   TOF  = %g [min]', TOF_constH_sub/60);
fprintf('\n Cruise Climb Time of Flight:               TOF  = %g [min]', TOF_CC_sub/60);
fprintf('\n\n -------------------------------------------------------------------- \n');
fprintf('\n Supersonic cruise: \n');
fprintf('\n Cruise Mach number: M  = %g', M_super_R);
fprintf('\n Cruise Velocity:    V  = %g [m/s] = %g [ft/s]', V_super_R, V_super_R*3.2808399);
fprintf('\n Cruise Altitude:    h  = %g [m] = %g [ft] \n', h_super_R, h_super_R*3.2808399);
fprintf('\n Constant Altitude Range:    R  = %g [km]   = %g [miles]', R_constH_super/1000, R_constH_super*0.000621371);
fprintf('\n Cruise Climb Range:         R  = %g [km]   = %g [miles] \n', R_CC_super/1000, R_CC_super*0.000621371);
fprintf('\n Constant Altitude Cruise Time of Flight:   TOF  = %g [min]', TOF_constH_super/60);
fprintf('\n Cruise Climb Time of Flight:               TOF  = %g [min]', TOF_CC_super/60);
fprintf('\n\n ============================================================================== \n');
%% ========================================================================
% Descent:
%--------------------------------------------------------------------------
V_descend = V_sub_R;                                  % [m/s] descent velocity
gamma_descend = -atand(1000/(3*5280));                % 3:1 rule
ROD = V_descend*sind(gamma_descend);                  % [m/s] rate of descent
%--------------------------------------------------------------------------
S_descend_sub = h_sub_R/abs(tand(gamma_descend));     % [m] range covered during descent
dt_descend_sub = S_descend_sub/V_descend;             % [s] time to descend
%--------------------------------------------------------------------------
S_descend_super = h_super_R/abs(tand(gamma_descend)); % [m] range covered during descent
dt_descend_super = S_descend_super/V_descend;         % [s] time to descend
%--------------------------------------------------------------------------
fprintf('\n\n ============================== Descent Results  ============================== \n');
fprintf('\n Descent from Subsonic Cruising Altitude: (h = %g [m])', h_sub_R);
fprintf('\n Descent Angle:         gamma = %g [deg] ', gamma_descend);
fprintf('\n Rate of Descent:         ROD = %g [m/s] = %g [ft/s]', ROD, ROD*3.2808399);
fprintf('\n Descent Velocity:          V = %g [m/s] = %g [ft/s]', V_descend, V_descend*3.2808399);
fprintf('\n\n -------------------------------------------------------------------- ');
fprintf('\n Range Covered During Descent:   R = %g [km] = %g [miles]', S_descend_sub/1000, S_descend_sub*0.000621371);
fprintf('\n Time to Descend:               dt = %g [min]', dt_descend_sub/60);
fprintf('\n\n -------------------------------------------------------------------- ');
fprintf('\n --------------------------------------------------------------------\n ');
fprintf('\n Descent from Supersonic Cruising Altitude: (h = %g [m])', h_super_R);
fprintf('\n Descent Angle:         gamma = %g [deg] ', gamma_descend);
fprintf('\n Rate of Descent:         ROD = %g [m/s] = %g [ft/s]', ROD, ROD*3.2808399);
fprintf('\n Descent Velocity:          V = %g [m/s] = %g [ft/s]', V_descend, V_descend*3.2808399);
fprintf('\n\n -------------------------------------------------------------------- ');
fprintf('\n Range Covered During Descent:   R = %g [km] = %g [miles]', S_descend_super/1000, S_descend_super*0.000621371);
fprintf('\n Time to Descend:               dt = %g [min]', dt_descend_super/60);
fprintf('\n\n ============================================================================== \n');
%% ========================================================================
% Landing:
%--------------------------------------------------------------------------
[~, ~, ~, rho_land, son_land, ~, ~, DELTA_land, THETA_land, ~] = ATMO(0, 'M');
CD_Land = Drag_Buildup_FUNKY(CL_Land, 0.3);
gamma_land = 3;                               % [deg] approach angle
V_stall = sqrt(2*W_land/(rho_TO*S_w*CL_max)); % [m/s] stall speed
V_approach = 1.2*V_stall;                     % [m/s] approach velocity
V_TD = 1.15*V_stall;                          % [m/s] touch down velocity

T_land = 0.01*T_max; % [N] landing thrust = idle thrust
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
mu = 0.5; % coefficient of friction

S_break = (W_land/(g*rho_land*S_w*(CD_Land-mu*CL_Land)))*log(1 + ((CD_Land-mu*CL_Land)*rho_land*S_w*V_TD^2)/(2*(-T_land + mu*W_land)));
%--------------------------------------------------------------------------
% total landing length:
S_land = S_approach + S_flare + S_Freeroll + S_break; % [m]
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
fprintf('\n Breaking Distance:               S_Breaking  = %g [m] = %g [ft]', S_break, S_break*3.2808399);
fprintf('\n\n Total Landing Distance:        S_Land = %g [m] = %g [ft]', S_land, S_land*3.2808399);
fprintf('\n\n FAR Landing Distance:          S_FAR  = %g [m] = %g [ft]', FAR_land, FAR_land*3.2808399);
fprintf('\n\n =============================================================================== \n');
%% ========================================================================
% total range and time of flight:
%--------------------------------------------------------------------------
% subsonic:
R_total_sub_1 = S_climb_sub + R_constH_sub + S_descend_sub; % [m]
R_total_sub_2 = S_climb_sub + R_CC_sub + S_descend_sub; % [m]

dt_total_sub_1 = dt_climb_sub + TOF_constH_sub + dt_descend_sub; % [s]
dt_total_sub_2 = dt_climb_sub + TOF_CC_sub + dt_descend_sub; % [s]
%--------------------------------------------------------------------------
% supersonic:
R_total_super_1 = S_climb_super + R_constH_super + S_descend_super; % [m]
R_total_super_2 = S_climb_super + R_CC_super + S_descend_super; % [m]

dt_total_super_1 = dt_climb_super + TOF_constH_super + dt_descend_super; % [s]
dt_total_super_2 = dt_climb_super + TOF_CC_super + dt_descend_super; % [s]
%--------------------------------------------------------------------------
fprintf('\n\n ============================== Total Range and Time of Flight ============================== \n');
fprintf('\n Total Range (including climb and descent) for Subsonic Constant Altitude Cruise:');
fprintf('\n -------------------------------------------------------------------------------- ');
fprintf('\n Cruise Altitude: h  = %g [m] = %g [ft] \n', h_sub_R, h_sub_R*3.2808399);
fprintf('\n Range Covered During Climb:        R_climb   = %g [km] = %g [miles]', S_climb_sub/1000, S_climb_sub*0.000621371);
fprintf('\n Constant Altitude Cruise Range:    R_cruise  = %g [km] = %g [miles]', R_constH_sub/1000, R_constH_sub*0.000621371);
fprintf('\n Range Covered During Descent:      R_descend = %g [km] = %g [miles]', S_descend_sub/1000, S_descend_sub*0.000621371);
fprintf('\n\n Total Range:  R = %g [km] = %g [miles] \n', R_total_sub_1/1000, R_total_sub_1*0.000621371);
fprintf('\n Time to Climb:    dt_climb   = %g [min]', dt_climb_sub/60);
fprintf('\n Time to Cruise:   dt_cruise  = %g [min]', TOF_constH_sub/60);
fprintf('\n Time to Descend:  dt_descend = %g [min]', dt_descend_sub/60);
fprintf('\n\n Total Time of Flight:   dt = %g [min]', dt_total_sub_1/60);
fprintf('\n                         dt = %g [hrs]', dt_total_sub_1/3600);
fprintf('\n\n ===================================================================================== \n');
fprintf('\n Total Range (including climb and descent) for Subsonic Cruise Climb:');
fprintf('\n -------------------------------------------------------------------- ');
fprintf('\n Cruise Altitude: h  = %g [m] = %g [ft] \n', h_sub_R, h_sub_R*3.2808399);
fprintf('\n Range Covered During Climb:        R_climb   = %g [km] = %g [miles]', S_climb_sub/1000, S_climb_sub*0.000621371);
fprintf('\n Cruise Climb Range:                R_cruise  = %g [km] = %g [miles]', R_CC_sub/1000, R_CC_sub*0.000621371);
fprintf('\n Range Covered During Descent:      R_descend = %g [km] = %g [miles]', S_descend_sub/1000, S_descend_sub*0.000621371);
fprintf('\n\n Total Range:  R = %g [km] = %g [miles] \n', R_total_sub_2/1000, R_total_sub_2*0.000621371);
fprintf('\n Time to Climb:    dt_climb   = %g [min]', dt_climb_sub/60);
fprintf('\n Time to Cruise:   dt_cruise  = %g [min]', TOF_CC_sub/60);
fprintf('\n Time to Descend:  dt_descend = %g [min]', dt_descend_sub/60);
fprintf('\n\n Total Time of Flight:   dt = %g [min]', dt_total_sub_2/60);
fprintf('\n                         dt = %g [hrs]', dt_total_sub_2/3600);
fprintf('\n\n ===================================================================================== ');
fprintf('\n ===================================================================================== \n');
fprintf('\n Total Range (including climb and descent) for Supersonic Constant Altitude Cruise:');
fprintf('\n ---------------------------------------------------------------------------------- ');
fprintf('\n Cruise Altitude: h  = %g [m] = %g [ft] \n', h_super_R, h_super_R*3.2808399);
fprintf('\n Range Covered During Climb:        R_climb   = %g [km] = %g [miles]', S_climb_super/1000, S_climb_super*0.000621371);
fprintf('\n Constant Altitude Cruise Range:    R_cruise  = %g [km] = %g [miles]', R_constH_super/1000, R_constH_super*0.000621371);
fprintf('\n Range Covered During Descent:      R_descend = %g [km] = %g [miles]', S_descend_super/1000, S_descend_super*0.000621371);
fprintf('\n\n Total Range:  R = %g [km] = %g [miles] \n', R_total_super_1/1000, R_total_super_1*0.000621371);
fprintf('\n Time to Climb:    dt_climb   = %g [min]', dt_climb_super/60);
fprintf('\n Time to Cruise:   dt_cruise  = %g [min]', TOF_constH_super/60);
fprintf('\n Time to Descend:  dt_descend = %g [min]', dt_descend_super/60);
fprintf('\n\n Total Time of Flight:   dt = %g [min]', dt_total_super_1/60);
fprintf('\n                         dt = %g [hrs]', dt_total_super_1/3600);
fprintf('\n\n ===================================================================================== \n');
fprintf('\n Total Range (including climb and descent) for Supersonic Cruise Climb:');
fprintf('\n ---------------------------------------------------------------------- ');
fprintf('\n Cruise Altitude: h  = %g [m] = %g [ft] \n', h_super_R, h_super_R*3.2808399);
fprintf('\n Range Covered During Climb:        R_climb   = %g [km] = %g [miles]', S_climb_super/1000, S_climb_super*0.000621371);
fprintf('\n Cruise Climb Range:                R_cruise  = %g [km] = %g [miles]', R_CC_super/1000, R_CC_super*0.000621371);
fprintf('\n Range Covered During Descent:      R_descend = %g [km] = %g [miles]', S_descend_super/1000, S_descend_super*0.000621371);
fprintf('\n\n Total Range:  R = %g [km] = %g [miles] \n', R_total_super_2/1000, R_total_super_2*0.000621371);
fprintf('\n Time to Climb:    dt_climb   = %g [min]', dt_climb_super/60);
fprintf('\n Time to Cruise:   dt_cruise  = %g [min]', TOF_CC_super/60);
fprintf('\n Time to Descend:  dt_descend = %g [min]', dt_descend_super/60);
fprintf('\n\n Total Time of Flight:   dt = %g [min]', dt_total_super_2/60);
fprintf('\n                         dt = %g [hrs]', dt_total_super_2/3600);
fprintf('\n\n ============================================================================================ \n');

%% ========================================================================
% Operational Envelope
%--------------------------------------------------------------------------
%WEIGHT  = 112200*4.44822;   % [N] operating weight
q_max = 86184.414;          % [Pa] max. dyn. pressure (p.104 Nicolai) = 1800 psf 
%--------------------------------------------------------------------------
% Altitude and Velocity range:
h = 0:500:25000; % [m]
V = 0:5:800;     % [m/s]
%--------------------------------------------------------------------------
for n = 1:length(h)

   [alt, press, temp, rho, a, ~, ~, ~, ~, ~] = ATMO(h(n), 'M');
 
for m = 1:length(V)
     M(m) = V(m)/a;
%--------------------------------------------------------------------------
% weight ratios (from Ben's weight script):
[weights] = AVD_Weight_Buildup_simple(M(m));
W_cruise_start = weights.W_accel_end;    % [lbf] weight at start of cruise
W_cruise_end   = weights.W_cruise_end;   % [lbf] weight at end of cruise
W_cruise_start = W_cruise_start*4.44822; % [N]
W_cruise_end   = W_cruise_end*4.44822;   % [N]
WEIGHT   = (W_cruise_start + W_cruise_end)/2; % [N] average cruise weight
%--------------------------------------------------------------------------
% Thrust available:
T_A_data = csvread('Thrust_NEW.csv'); % [lbf] per engine
T_A_data = T_A_data*4.44822*3;        % [N] all three engines
MM = 0:0.1:2;
hh = 0:1000:100000; % [ft]
hh = hh*0.3048;     % [m]

Thrust_A = interp2(MM,hh,T_A_data, M(m), h(n)); % [N] thrust available
%--------------------------------------------------------------------------
CL(n,m) = WEIGHT/(0.5*rho*V(m)^2*S_w); % steady level flight
%--------------------------------------------------------------------------
CD = Drag_Buildup_FUNKY(CL(n,m), M(m));
D = CD*0.5*rho*V(m)^2*S_w; % [N] Drag (sluf)

Ps(n,m) = V(m)*(Thrust_A - D)/WEIGHT; % Specific Excess Power
%--------------------------------------------------------------------------
end;
 V_stall_op(n) = sqrt(2*WEIGHT/(CL_max*rho*S_w)); %stall boundary
 M_stall_op(n) = V_stall_op(n)/a;
 V_q_lim(n) = sqrt(2*q_max/rho); % dynamic pressure limit
 M_q_lim(n) = V_q_lim(n)/a;
end;
%--------------------------------------------------------------------------
figure('Name','Operational Envelope','NumberTitle','off');
hold on;
vals = [1,1];
[C1, h1] = contour(M,h,Ps,vals);
set(h1, 'LineWidth', 2.5)
set(h1, 'Color', 'k')
title('Operational Envelope','FontSize',18);
xlabel('Mach Number','FontSize',12);
ylabel('Altitude (m)','FontSize',12);
xlim([0, 2.5]);
ylim([0, 2e4]);
plot (M_stall_op,h,':k','LineWidth',2.5);
plot (M_q_lim,h,'--k','LineWidth',2.5);
plot (M_sub_R,h_sub_R,'k+','LineWidth',2.5,'LineStyle','None','MarkerSize',16);
plot (M_super_R,h_super_R,'k*','LineWidth',2.5,'LineStyle','None','MarkerSize',16);
legend('Thrust Limit','Stall Boundary', 'Dyn. Press. Limit', 'Subsonic Cruise', 'Supersonic Cruise','Location','NorthWest');
%--------------------------------------------------------------------------
%% ========================================================================
% thrust required and thrust available:
%--------------------------------------------------------------------------
h = [h_sub_R, h_super_R]; % [m] Altitude 
M = 0.1:0.01:1.4;            % Mach #
%--------------------------------------------------------------------------
% thrust required:
for n = 1:length(h)
    
    [~, ~, ~, rho, a, ~, ~, ~, ~, ~ ] = ATMO(h(n), 'M');
    
    for m = 1:length(M)

    V = M(m)*a; % [m/s] velocity
%--------------------------------------------------------------------------
    % weight ratios (from Ben's weight script):
    [weights] = AVD_Weight_Buildup_simple(M(m));
    W_cruise_start = weights.W_accel_end;           % [lbf] weight at start of cruise
    W_cruise_end   = weights.W_cruise_end;          % [lbf] weight at end of cruise
    W_cruise_start = W_cruise_start*4.44822;        % [N]
    W_cruise_end   = W_cruise_end*4.44822;          % [N]
    WEIGHT   = (W_cruise_start + W_cruise_end)/2;   % [N] average cruise weight
%--------------------------------------------------------------------------
 	CL = WEIGHT/(0.5*rho*V^2*S_w);                  % steady level flight
    CD = Drag_Buildup_FUNKY(CL, M(m));              % Drag coefficient
    D = CD*0.5*rho*V^2*S_w;                         % [N] Drag (sluf)
    T_R(n,m) = D;                                   % [N] Thrust required
    T_A(n,m) = interp2(MM,hh,T_A_data, M(m), h(n)); % [N] thrust available

    end
end
%--------------------------------------------------------------------------
% at subsonic altitude:
figure('Name','T_r and T_a Subsonic Alt','NumberTitle','off');
hold on;
plot (M, T_R(1,:),'k','LineWidth',2.5);
plot (M, T_A(1,:),'--k','LineWidth',2.5);
title('Thrust Required and Available at subsonic cruising altitude','FontSize',18);
xlabel('Mach Number','FontSize',12);
ylabel('Thrust (N)','FontSize',12);
xlim([0.6, 1.4])
legend('Thrust Required','Thrust Available','Location','NorthEast');
%--------------------------------------------------------------------------
% at supersonic altitude:
figure('Name','T_r and T_a Supersonic Alt','NumberTitle','off');
hold on;
plot (M, T_R(2,:),'k','LineWidth',2.5);
plot (M, T_A(2,:),'--k','LineWidth',2.5);
title('Thrust Required and Available at supersonic cruising altitude','FontSize',18);
xlabel('Mach Number','FontSize',12);
ylabel('Thrust (N)','FontSize',12);
xlim([0.6, 1.4])
legend('Thrust Required','Thrust Available','Location','NorthEast');
%--------------------------------------------------------------------------