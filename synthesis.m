%{

    Integrated synthesis script for MAE 4351
    Aerospace senior design
%--------------------------------------------------------------------------
Team: ARROW
Team members: 
Shawn McCullough, Ben Holden, Nick Schulte, Rustin Farris, Christian Allen
%--------------------------------------------------------------------------
Last modified: 03/20/2018
% =========================================================================
%}
clear; clc; close all;
%% ========================================================================
% Mission inputs
%--------------------------------------------------------------------------
supersonic = 0; % Enter 1 if analyzing supersonic cruise, 0 for subsonic
%--------------------------------------------------------------------------
M_cr_sub = 0.9;              % Subsonic cruise Mach number
M_cr_super = 1.4;            % Supersonic cruise Mach number
range_sub = 9800*10^3;       % Subsonic range, (m)
range_super = 8800*10^3;     % Supersonic range, (m)
alt_sub_cr = 13000;          % Subsonic cruise altitude (m)
alt_super_cr = 15550;        % Supersonic cruise altitude (m)

M_max = 1.5; % From AS2 mission (need to update)
%--------------------------------------------------------------------------
cruise_altitude = round(alt_super_cr*3.2808399); % ft
altitudes = [0, 30000, cruise_altitude];         % ft (vor v-n diagram)
%--------------------------------------------------------------------------
if supersonic == 0
    alt_cr = alt_sub_cr;
    M_cr = M_cr_sub;
    R_cr = range_sub;
else 
    alt_cr = alt_super_cr;
    M_cr = M_cr_super;
    R_cr = range_super;
end
%% ========================================================================
% Empirical inputs
%--------------------------------------------------------------------------
TSFC = 1.0;    % Empirical Placeholder for SFC based on Sadraey Table 4.6 for turbojet
LD_cruise = 9; % Cruise lift/drag from Fig 5.3 Nicolai
CL_max = 2.6;  % Placeholder, max CL

%% ========================================================================
% Interdisciplinary inputs (Design inputs) 
%--------------------------------------------------------------------------
AR = 10;            % Placeholder, unswept aspect ratio
AR_low = AR;        % Low speed, unswept AR
TR = 0.44;          % Placeholder wing taper ratio
tmax = 2.3;         % Maximum thickness, based on AS2 cabin dimensions (m)
e = 0.85;           % Oswald
K = 1/(pi*AR*e);	% Drag K factor
ne = 4;             % number of engines
SM = 1;             % static margin

%% ========================================================================
% Solution Space
%--------------------------------------------------------------------------
[design_point] = Solution_Space_OFW_integrated(AR);

WingLoading   = design_point(1); % W/S (lbf/ft^2)
ThrustLoading = design_point(2); % T/W (lbf/lbf)

%% ========================================================================
% Weights
%--------------------------------------------------------------------------
[~, ~, ~, ~, a] = ATMO(cruise_altitude, 'M');
num_pass = 19;                  % number of passengers
num_crew = 4;                   % number of crew members
R_miles = R_cr*0.000621;        % range in miles
V_cr = M_cr*a;                  % cruise velocity, (m/s) 
V_cr_mph = V_cr*2.23694;        % Cruise vel (mph) for weight script
%--------------------------------------------------------------------------
[weights, wt_frac] = Weight_Buildup(num_pass, num_crew, V_cr_mph, M_cr, R_miles, TSFC, LD_cruise);

MTOW = 4.44822*weights.W_to;    % Max takeoff weight, (N)
W_to_end = MTOW*wt_frac.WF_to;  % Wt at end of TO, start of climb (N)

W_climb_end = MTOW*wt_frac.WF_to*wt_frac.WF_climb*wt_frac.WF_accel; % Wt at end of climb, start of cruise (N)

W_climb_avg = 0.5*(W_to_end + W_climb_end);

W_cruise_start = 4.44822*W_climb_end; % Wt at beginning of cruise (N)
W_cruise_end = MTOW*wt_frac.WF_to*wt_frac.WF_climb*wt_frac.WF_accel*wt_frac.WF_cruise; % Wt at end of cruise (N)
W_cruise_avg = 0.5*(W_climb_end + W_cruise_end); % Average cruise wt (N)

W_descend_end = MTOW*wt_frac.WF_to*wt_frac.WF_climb*wt_frac.WF_accel*wt_frac.WF_cruise*wt_frac.WF_des; % Wt at end of descent (lbf)
W_descend_avg = 0.5*(W_cruise_end + W_descend_end); % Average descent wt (N)

W_land = 4.44822*weights.W_land; % Landing Weight (N)

%% ========================================================================
% required wing area and takeoff thrust
%--------------------------------------------------------------------------
S = weights.W_to/WingLoading; % (ft^2) for V-n diagram
Sref = S*0.092903;            % Wing area in (m^2) for everything else

T_max = 4*ThrustLoading*MTOW; % (N) required take off thrust

%% ========================================================================
% more Interdisciplinary inputs
%--------------------------------------------------------------------------
b_ft = sqrt(AR*S);  % (ft) for wing loading 
b = b_ft*0.3048;    % Wing span, (m) for everything else
%--------------------------------------------------------------------------
% Calculate sweep:
M_perp = 0.7;                       % perp Mach #, variable to iterate
sweep_deg = acosd(M_perp/M_cr);     % This has to be limited to 70 ish degrees!
if sweep_deg > 70                   % Limit the sweep angle
    sweep_deg = 70;
end
sweep_rad = sweep_deg*pi/180;       % Convert to radians for formulas
b_swept = b*cosd(sweep_deg);        % Span at sweep angle
AR_swept = b_swept^2/AR;            % Swept aspect ratio
%--------------------------------------------------------------------------
if M_cr > 0.69
    AR = AR_swept;
    b = b_swept;
end
%--------------------------------------------------------------------------

%% ========================================================================
% display initial design parameters:
%--------------------------------------------------------------------------
fprintf('\n\n ============================= Initial Design Parameters ============================= \n');
fprintf('\n Required wing area: S  = %g [m^2] = %g [ft^2] ', Sref, S);
fprintf('\n Wing span:          b   = %g [m] = %g [ft]', b, b_ft);
fprintf('\n Cruise Wing sweep:  Lambda = %g [deg] \n', sweep_deg);
fprintf('\n -------------------------------------------------------------------------------- ');
fprintf('\n Required takeoff thrust:   T  = %g [N] = %g [lbf] ', T_max, T_max*0.22480902);
fprintf('\n Max takeoff weight:        MTOW  = %g [N] = %g [lbf] ', MTOW, MTOW*0.22480902);
fprintf('\n\n ===================================================================================== \n');

%% ========================================================================
% Preliminary Analysis ????????????????????????????????????????????????????
%--------------------------------------------------------------------------
% Structures, W&B

%% ========================================================================
% V-n diagram
%--------------------------------------------------------------------------
Vn_Diagram(convforce(MTOW,'N','lbf'), Sref, altitudes, M_cr, M_max, CL_max);
[max_load, min_load] = Wing_Loading(b_ft, MTOW, TR); % MTOW is already in Newtons!
%[max_load, min_load] = Wing_Loading(b_ft, convforce(MTOW, 'lbf', 'N'), TR);

%% ========================================================================
% Takeoff Phase
%--------------------------------------------------------------------------
alt_TO = 0; % Airport altitude [m]
[~, ~, ~, rho_TO, son_TO, ~, ~, ~, ~, ~] = ATMO(alt_TO, 'M');
V_TO = 86;          % Placeholder, [m/s]
M_TO = V_TO/son_TO;	% Mach number for takeoff phase
%--------------------------------------------------------------------------
% aero:
CL_TO = 2*MTOW/(0.5*rho_TO*V_TO^2*Sref);   % N/N
[CD_TO, CD0_TO] = aerofunk_drag_2(alt_TO, M_TO, Sref, CL_TO, SM, AR, TR);

%--------------------------------------------------------------------------
% performance:
[V_stall, V_TO, S_G, S_TO, BFL] = perf_takeoff(ne, V_TO, T_max, alt_TO, MTOW, Sref, CL_max, CD_TO, K);

%--------------------------------------------------------------------------
% S&C:
[CMa, Cl_beta, Cn_beta, Cn_dr, S_VT, l_VT] = stability(M_TO, AR, 24, Sref, b, TR, CL_TO, SM, 'Takeoff');

%% ========================================================================
% Climb Phase
%--------------------------------------------------------------------------
CL_climb = CL_max/(1.25^2);           % climb lift coefficient
V_sub_R = M_cr_sub * 295.07;          % Cruise velocity in m/s
V = linspace(V_stall,V_sub_R)';       % [m/s] velocity vector
V_climb = V_sub_R*0.85;               % [m/s] climb velocity
M_climb = V_climb/son_TO;             % climb Mach number
%[T, T_max] = Thrust(alt_cr, M_cr);
[~, ~, ~, ~, a] = ATMO(alt_cr, 'm');
V_sub = M_cr_sub*a;                   % cruise velocity, m/s 

if M_cr > 1
    alt_climb = alt_sub_cr;              
else
    alt_climb = alt_super_cr;
end

AoA_climb = 0; % fix this!
[CD_climb, CD0_climb] = aerofunk_drag_2(alt_climb, M_climb, Sref, CL_climb, SM, AR, TR);

[ROC, gamma_climb, S_climb, dt_climb] = perf_climb(CL_climb, CD_climb, V_stall, V_sub, rho_TO, T_max, W_climb_avg, alt_climb, Sref);

%% ========================================================================
% Cruise Phase
%--------------------------------------------------------------------------
[~, ~, ~, rho_cr, son_cr, ~, ~, ~, ~, ~] = ATMO(alt_cr, 'M');
V_cr = M_cr*son_cr; % [m/s]

CL_cr = W_cruise_start/(Sref*0.5*rho_cr*V_cr^2);     % lift coefficient cruise

[R_constH, R_CC, TOF_constH, TOF_CC] = perf_cruise(M_cr, alt_cr, W_cruise_start, W_cruise_end, Sref, SM, AR, TSFC, TR);

%--------------------------------------------------------------------------
% S&C:
[CMa, Cl_beta, Cn_beta, Cn_dr, S_VT, l_VT] = stability(M_cr, AR, 24, Sref, b, TR, CL_cr, SM, 'Cruise');

%% ========================================================================
% Descent Phase
%--------------------------------------------------------------------------
V_descend = V_cr;
[S_descend, dt_descend] = perf_descent(alt_cr, V_descend);

%% ========================================================================
% Landing Phase
%--------------------------------------------------------------------------
M_land = 0.3;             % landing Mach number
alt_land = 0;

[S_land, FAR_land] = perf_land(M_land, alt_land, Sref, AR_low, W_land, CL_max, T_max, TR, SM);

%--------------------------------------------------------------------------
% S&C:
[CMa, Cl_beta, Cn_beta, Cn_dr, S_VT, l_VT] = stability(M_land, AR, 24, Sref, b, TR, CL_max, SM, 'Landing');


%% ========================================================================
% Total Performance Summary:
%--------------------------------------------------------------------------
% total range and time of flight:
R_total_1 = S_climb + R_constH + S_descend; % [m]
R_total_2 = S_climb + R_CC + S_descend; % [m]

dt_total_1 = dt_climb + TOF_constH + dt_descend; % [s]
dt_total_2 = dt_climb + TOF_CC + dt_descend; % [s]
%--------------------------------------------------------------------------
fprintf('\n\n ============================== Total Range and Time of Flight ============================== \n');
fprintf('\n Total Range (including climb and descent) for Constant Altitude Cruise:');
fprintf('\n -------------------------------------------------------------------------------- ');
fprintf('\n Cruise Altitude: h  = %g [m] = %g [ft] \n', alt_cr, alt_cr*3.2808399);
fprintf('\n Range Covered During Climb:        R_climb   = %g [km] = %g [miles]', S_climb/1000, S_climb*0.000621371);
fprintf('\n Constant Altitude Cruise Range:    R_cruise  = %g [km] = %g [miles]', R_constH/1000, R_constH*0.000621371);
fprintf('\n Range Covered During Descent:      R_descend = %g [km] = %g [miles]', S_descend/1000, S_descend*0.000621371);
fprintf('\n\n Total Range:  R = %g [km] = %g [miles] \n', R_total_1/1000, R_total_1*0.000621371);
fprintf('\n Time to Climb:    dt_climb   = %g [min]', dt_climb/60);
fprintf('\n Time to Cruise:   dt_cruise  = %g [min]', TOF_constH/60);
fprintf('\n Time to Descend:  dt_descend = %g [min]', dt_descend/60);
fprintf('\n\n Total Time of Flight:   dt = %g [min]', dt_total_1/60);
fprintf('\n                         dt = %g [hrs]', dt_total_1/3600);
fprintf('\n\n ===================================================================================== \n');
%% ========================================================================
% operational envelope:
%--------------------------------------------------------------------------
op_envelope(W_cruise_avg, T_max, Sref, SM, AR, TR, CL_max)

%% ========================================================================
% Cost analysis:
%--------------------------------------------------------------------------
W_A = 0; % what is this?
[ RTDE_Cost ] = costfunky( weights.W_empty, V_cr, T_max, M_max, MTOW, weights.W_fuel, R_total_1*0.000621371, dt_climb/3600, dt_descend/3600, TOF_constH/3600, W_A, ne);

%% >>>>>>> .r145 <----what is this??
