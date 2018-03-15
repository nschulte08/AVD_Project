clear; clc; close all;

%% Mission inputs
supersonic = 1;         % Enter 1 if analyzing supersonic cruise, 0 for sub
cruise_altitude = 51000; % ft
altitudes = [0, 30000, cruise_altitude]; % ft
alt_sub_cr = 13000;             % Subsonic cruise altitude (m)
alt_super_cr = 15550;           % Supersonic cruise altitude (m)
M_cr_sub = 0.95;
M_cr_super = 1.4;
range_sub = 9800*10^3;       % Subsonic range, (m)
range_super = 8800*10^3;     % Supersonic range, (m)
if supersonic == 0
    alt_cr = alt_sub_cr;
    M_cr = M_cr_sub;
    R_cr = range_sub;
else 
    alt_cr = alt_super_cr;
    M_cr = M_cr_super;
    R_cr = range_super;
end  
num_pass = 19;                  % number of passengers
num_crew = 4;                   % number of crew members
R_miles = R_cr*.000621;         % range in miles
[~, ~, ~, ~, a] = ATMO(cruise_altitude, 'M');
V_cr = M_cr * a;                % cruise velocity, (m/s) 
V_cr_mph = V_cr*2.23694;            % Cruise vel (mph) for weight script


%% Empirical inputs
SFC = 1.0;    % Empirical Placeholder for SFC based on Sadraey Table 4.6 for turbojet
LD_cruise = 9; % Cruise lift/drag from Fig 5.3 Nicolai

%% Interdisciplinary (Design) inputs
S = 15700;          % Placeholder wing area, (ft^2) for V-n diagram
Sref = S * 0.092903; % Wing area in (m^2) for everything else
b_ft = 120;            % Placeholder wingspan, (ft) for wing loading 
b = b_ft * 0.3048;  % Wing span, (m) for everything else
tmax = 2.3;         % Maximum thickness, based on AS2 cabin dimensions (m)
TR = 0.44;      % Placeholder wing taper ratio
M_max = 1.5;        % From AS2 mission
CL_max = 2.6;       % Placeholder, max CL
AR = 10;            % Placeholder, unswept aspect ratio
AR_low = AR;        % Low speed, unswept AR

[~, T_max] = Thrust(alt_cr, M_cr);


% Calculate sweep:
M_perp = 0.7;                       % perp Mach #, variable to iterate
sweep_deg = acosd(M_perp/M_cr);        % This has to be limited to 70 ish degrees!
if sweep_deg > 70                   % Limit the sweep angle
    sweep_deg = 70;
end
sweep_rad = sweep_deg*pi/180;       % Convert to radians for formulas

b_swept = b*cosd(sweep_deg);        % Span at sweep angle
AR_swept = b_swept^2 / AR;          % Swept aspect ratio

if M_cr > .69
    AR = AR_swept;
    b = b_swept;
end

e = 0.85;           % Oswald
Kfact = 1/(pi*AR*e);	% Drag K factor


%% Preliminary Analysis
% Structures, W&B

%% Phase-independent function calls

[weights, wt_frac] = Weight_Buildup(num_pass, num_crew, V_cr_mph, M_cr, R_miles, SFC, LD_cruise);

MTOW = 4.44822 * weights.W_to;    % Max takeoff weight, converted to N
W_to_end = MTOW * wt_frac.WF_to;                                            % Wt at end of TO, start of climb (lbf)

W_climb_end = MTOW * wt_frac.WF_to * wt_frac.WF_climb * wt_frac.WF_accel;     % Wt at end of climb, start of cruise (lbf)
W_climb_avg = 0.5 * (W_to_end + W_climb_end);

W_cruise_end = MTOW * wt_frac.WF_to * wt_frac.WF_climb * ...
    wt_frac.WF_accel * wt_frac.WF_cruise;                       % Wt at end of cruise
W_cruise_avg = 0.5 * (W_climb_end+W_cruise_end);                % Average cruise wt (lbf)

W_descend_end = MTOW * wt_frac.WF_to * wt_frac.WF_climb * ...
    wt_frac.WF_accel * wt_frac.WF_cruise * wt_frac.WF_des;      %Wt at end of descent (lbf)
W_descend_avg = 0.5 * (W_cruise_end + W_descend_end);

W_land = weights.W_land;                                        % Landing Weight (lbf)

Vn_Diagram(MTOW, Sref, altitudes, M_cr, M_max, CL_max);
[max_load, min_load] = Wing_Loading(b_ft, convforce(MTOW, 'lbf', 'N'), TR);

% Flight Phase Analysis

%% Takeoff & Climb Phase
alt_TO = 0;                 % Airport altitude
[~, ~, ~, rho_TO, son_TO, ~, ~, ~, ~, SIGMA_TO] = ATMO(alt_TO, 'M');
V_TO = 86;                                           % Placeholder, [m/s]
M_TO = V_TO/son_TO;                                  % Mach number for takeoff phase


% Calculate lift coefficient for takeoff as a function of MTOW
CL_TO = 2 * W_climb_avg / (0.5 * rho_TO * V_TO^2 * Sref);   % N / N

alpha_TO = 1.0;             % Placeholder, need to find from CLa and CL from aerofunk_lift

% Calculate CD and CD0 for takeoff as a function of CL and M
[CD_TO, CD0_TO] = aerofunk_drag(alt_TO, M_TO, Sref, CL_TO, alpha_TO, AR, TR);
[V_stall, V_TO, S_G, S_TO, BFL] = perf_takeoff(T_max, alt_TO, MTOW, Sref, CL_max, CD0_TO, Kfact);


% Put S&C stuff here


%% Climb Phase
CL_climb = CL_max/(1.25^2);         % climb lift coefficient
V_sub_R = M_cr_sub * 295.07;        % Cruise velocity in m/s
V = linspace(V_stall,V_sub_R)';                    % [m/s] velocity vector
V_climb = floor(V(ind));                           % [m/s] climb velocity
M_climb = V_climb/son_TO;                          % climb Mach number
[T, T_max] = Thrust(alt_cr, M_cr);
[~, ~, ~, ~, a] = ATMO(alt_cr_sub, 'm');
V_sub = M_cr_sub * a;                % cruise velocity, m/s 

if M_cr > 1
    alt_climb = alt_sub_cr;              
else
    alt_climb = alt_super_cr;
end

[CD_climb, CD0_climb] = aerofunk_drag(alt_climb, M_climb, Sref, CL_climb, alpha, AR, TR);

[ROC, gamma_climb, S_climb, dt_climb] = perf_climb(CD_climb, V_stall, V_sub, rho_TO, T_max, W_climb_avg, alt_climb, V_climb);

%% Cruise Phase
W_cr_start = W_climb_end;           % Wt at beginning of cruise
W_cr_end = W_cruise_end;            % Wt at end of cruise
    
alpha_cr = 1.0;     %placeholder, need to calculate alpha reqd based on CLa and cruise CL

%Calculate subsonic, steady level LD max
LD_max = 0.5*sqrt((pi*AR*e)/CD0_sub);                 

[R_constH, R_CC, TOF_constH, TOF_CC] = perf_cruise(M_cr, alt_cr, W_cr_start, W_cr_end, Sref, alpha_cr, AR, TSFC);


% Put S&C stuff here



%% Descent Phase

[S_descend, dt_descend] = perf_descent(alt_sub_cr, M_cr);


%% Landing Phase
alt_land = alt_TO;

[S_land, FAR_land] = perf_land(alt_land, Sref, W_descend_avg, AR_low, W_land, CL_max, T_max);


% Put S&C stuff here


%% Total Performance Summary:
%% ========================================================================
% total range and time of flight:
%--------------------------------------------------------------------------
R_total_1 = S_climb + R_constH + S_descend; % [m]
R_total_2 = S_climb + R_CC + S_descend; % [m]

dt_total_1 = dt_climb + TOF_constH + dt_descend; % [s]
dt_total_2 = dt_climb + TOF_CC + dt_descend; % [s]

%--------------------------------------------------------------------------
fprintf('\n\n ============================== Total Range and Time of Flight ============================== \n');
fprintf('\n Total Range (including climb and descent) for Constant Altitude Cruise:');
fprintf('\n -------------------------------------------------------------------------------- ');
fprintf('\n Cruise Altitude: h  = %g [m] = %g [ft] \n', alt_cr, alt_cr*3.2808399);
fprintf('\n Range Covered During Climb:        R_climb   = %g [km] = %g [miles]', S_climb/1000, S_climb_sub*0.000621371);
fprintf('\n Constant Altitude Cruise Range:    R_cruise  = %g [km] = %g [miles]', R_constH/1000, R_constH_sub*0.000621371);
fprintf('\n Range Covered During Descent:      R_descend = %g [km] = %g [miles]', S_descend/1000, S_descend_sub*0.000621371);
fprintf('\n\n Total Range:  R = %g [km] = %g [miles] \n', R_total_1/1000, R_total_1*0.000621371);
fprintf('\n Time to Climb:    dt_climb   = %g [min]', dt_climb/60);
fprintf('\n Time to Cruise:   dt_cruise  = %g [min]', TOF_constH/60);
fprintf('\n Time to Descend:  dt_descend = %g [min]', dt_descend/60);
fprintf('\n\n Total Time of Flight:   dt = %g [min]', dt_total_1/60);
fprintf('\n                         dt = %g [hrs]', dt_total_1/3600);
fprintf('\n\n ===================================================================================== \n');




%=======
%% Takeoff
% Aerodynamics function calls
% Propulsion function calls
% Performance function calls
% Structures, W&B function calls
% Stability & control function calls

%% Climb
% Performance function calls

%% Cruise
% Aerodynamics function calls
% Propulsion function calls
% Performance function calls
% Stability & control function calls

%% Descent
% Performance function calls

%% Landing
% Aerodynamics function calls
% Propulsion function calls
% Performance function calls
% Structures, W&B function calls
% Stability & control function calls
%>>>>>>> .r145

