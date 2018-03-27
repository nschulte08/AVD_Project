%{
    Integrated synthesis script for MAE 4351
    Aerospace senior design
%--------------------------------------------------------------------------
Team: ARROW
Team members: 
Shawn McCullough, Ben Holden, Nick Schulte, Rustin Farris, Christian Allen
%--------------------------------------------------------------------------
Last modified: 03/26/2018
% =========================================================================
%}
clear; clc; close all;
%% ========================================================================
% Mission inputs
%--------------------------------------------------------------------------
supersonic = 1; % Enter 1 if analyzing supersonic cruise, 0 for subsonic
%--------------------------------------------------------------------------
M_cr_sub = 0.9;              % Subsonic cruise Mach number
M_cr_super = 1.4;            % Supersonic cruise Mach number
range_sub = 9800e3;          % Subsonic range, (m)
range_super = 8800e3;        % Supersonic range, (m)
alt_sub_cr = 13000;          % Subsonic cruise altitude (m)
alt_super_cr = 15500;        % Supersonic cruise altitude (m)
M_max = 1.5;                 % From AS2 mission (need to update)
num_pass = 19;               % Number of passengers
num_crew = 4;                % Number of crew members
%--------------------------------------------------------------------------
alt_TO = 0; % takeoff Airport altitude [m]
alt_land = 0; % landing Airport altitude [m]
%--------------------------------------------------------------------------
if supersonic == 0
    alt_cr = alt_sub_cr; % [m]
    M_cr = M_cr_sub;
    R_cr = range_sub; % [m]
    [~,TSFC] = Propulsion(M_cr_sub,alt_sub_cr); % TSFC = [1/hr]
elseif supersonic == 1
    alt_cr = alt_super_cr; % [m]
    M_cr = M_cr_super;
    R_cr = range_super; % [m]
    [~,TSFC] = Propulsion(M_cr_super,alt_super_cr); % TSFC = [1/hr]
end
%% ========================================================================
% Empirical inputs
%--------------------------------------------------------------------------
%TSFC = 0.9;     % [1/hr] Empirical Placeholder for SFC (based on Sadraey Table 4.6 for turbojet = 1.0, turbofan = ???, )
LD_cruise = 9;  % Cruise lift/drag from Fig 5.3 Nicolai
CL_max = 2.9;   % Placeholder, max CL

%% ========================================================================
% Interdisciplinary inputs (Design inputs) 
%--------------------------------------------------------------------------
% wing geometry:
AR_unswept = 8;                         % Unswept aspect ratio
AR_lowspeed = AR_unswept;               % Low speed, unswept AR
TR = 0.40;                              % Wing taper ratio
tmax = 2.3;                             % Maximum thickness, based on AS2 cabin dimensions (m)
tcmax = 0.16;                           % T/c max; variable to iterate
c_r = tmax/tcmax;                       % [m] Root chord = max thickness / tcmax ratio
c_t = TR*c_r;                           % [m] Root chord = max thickness / tcmax ratio
b_unswept = (AR_unswept/2)*(c_r + c_t); % [m] wing span
Sref = b_unswept^2/AR_unswept;          % [m^2] wing area
%--------------------------------------------------------------------------
% aerodynamics and S&C: 
e = 0.85;           % Oswald
K = 1/(pi*AR_unswept*e);	% Drag K factor
SM = 1;             % static margin
%--------------------------------------------------------------------------
% propulsion:
ne = 4; % number of engines

%% ========================================================================
% Solution Space
%--------------------------------------------------------------------------
[design_point] = Solution_Space_OFW_integrated(AR_unswept, CL_max, e, alt_sub_cr, M_cr_sub, alt_super_cr, M_cr_super, ne);

WingLoading   = design_point(1); % W/S (lbf/ft^2)
ThrustLoading = design_point(2); % T/W (lbf/lbf)

MTOW = convforce(WingLoading*Sref*convlength(1,'m','ft')^2, 'lbf', 'N'); % Max takeoff weight, N
T_max_required = ThrustLoading*MTOW; % Required thrust for takeoff, N

%% ========================================================================
% Weights
%--------------------------------------------------------------------------
[~, ~, ~, ~, son_climb, ~, ~, ~, ~, ~] = ATMO(alt_super_cr, 'M');
V_cr = M_cr*son_climb; % cruise velocity, (m/s) 
%--------------------------------------------------------------------------
[weights, wt_frac] = Weight_Buildup(convforce(MTOW,'N','lbf'), num_pass, num_crew, convvel(V_cr,'m/s','mph'), M_cr, convlength(R_cr,'m','mi'), TSFC, LD_cruise);

W_to_end = MTOW*wt_frac.WF_to;  % Wt at end of TO, start of climb (N)

W_climb_end = MTOW*wt_frac.WF_to*wt_frac.WF_climb*wt_frac.WF_accel; % Wt at end of climb, start of cruise (N)
W_climb_avg = 0.5*(W_to_end + W_climb_end);

W_cruise_start = W_climb_end; % Wt at beginning of cruise (N)
W_cruise_end = MTOW*wt_frac.WF_to*wt_frac.WF_climb*wt_frac.WF_accel*wt_frac.WF_cruise; % Wt at end of cruise (N)
W_cruise_avg = 0.5*(W_climb_end + W_cruise_end); % Average cruise wt (N)

W_descend_end = MTOW*wt_frac.WF_to*wt_frac.WF_climb*wt_frac.WF_accel*wt_frac.WF_cruise*wt_frac.WF_des; % Wt at end of descent (lbf)
W_descend_avg = 0.5*(W_cruise_end + W_descend_end); % Average descent wt (N)

W_land = convforce(weights.W_land,'lbf','N'); % Landing Weight (N)

%% ========================================================================
% Calculate sweep for cruise:
M_perp = 0.7;                        % perp Mach #, variable to iterate(?)
sweep_deg = acosd(M_perp/M_cr);      % This has to be limited to 70 ish degrees!
if sweep_deg > 70                    % Limit the sweep angle
    sweep_deg = 70;
end
sweep_rad = sweep_deg*pi/180;        % Convert to radians for formulas
b_swept = b_unswept*cosd(sweep_deg); % Span at sweep angle [m]
AR_swept = b_swept^2/Sref;           % Swept aspect ratio
%--------------------------------------------------------------------------
% if M_cr > 0.69 % what is this for???
%     AR = AR_swept;
%     b_unswept = b_swept; % [m]
% end
%--------------------------------------------------------------------------

%% ========================================================================
% Display initial design parameters:
%--------------------------------------------------------------------------
fprintf('\n\n ============================= Initial Design Parameters ============================= \n');
fprintf('\n Required wing area:         S  = %g [m^2] = %g [ft^2] ', Sref, Sref*convlength(1,'m','ft')^2);
fprintf('\n Total unswept wing span:    b  = %g [m] = %g [ft]', b_unswept, convlength(b_unswept,'m','ft')/cosd(sweep_deg));
fprintf('\n Effective wing span:        b_eff  = %g [m] = %g [ft]', b_swept, convlength(b_swept,'m','ft'));
fprintf('\n Cruise wing sweep:          Lambda = %g [deg] \n', sweep_deg);
fprintf('\n -------------------------------------------------------------------------------- ');
fprintf('\n Required takeoff thrust:   T  = %g [N] = %g [lbf] ', T_max_required, convforce(T_max_required,'N','lbf'));
fprintf('\n Max takeoff weight:        MTOW  = %g [N] = %g [lbf] ', MTOW, convforce(MTOW,'N','lbf'));
fprintf('\n\n ===================================================================================== \n');

%% ========================================================================
% operational envelope:
%--------------------------------------------------------------------------
%op_envelope(W_cruise_avg, Sref, SM, b_unswept, TR, CL_max, ne)

%% ========================================================================
% V-n diagram and wing loading
%--------------------------------------------------------------------------
altitudes = [0, convlength(alt_sub_cr,'m','ft'), convlength(alt_super_cr,'m','ft')]; % array of key altitudes for V-n diagram (ft)
Vn_Diagram(convforce(MTOW,'N','lbf'), Sref*convlength(1,'m','ft')^2, altitudes, M_cr, M_max, CL_max);
[max_load, min_load] = Wing_Loading(b_unswept, MTOW, TR); 

%% ========================================================================
% Takeoff Phase
%--------------------------------------------------------------------------
[~, ~, ~, rho_TO, son_TO, ~, ~, ~, ~, ~] = ATMO(alt_TO, 'M');

V_stall = sqrt(2*MTOW/(rho_TO*Sref*CL_max)); % [m/s]
V_TO = 1.1*V_stall; % [m/s]
M_TO = V_TO/son_TO;	% Mach number for takeoff phase
%--------------------------------------------------------------------------
% aero:
CL_TO = 2*MTOW/(0.5*rho_TO*V_TO^2*Sref);   % N/N
[CD_TO, CD0_TO] = aerofunk_drag_2(alt_TO, M_TO, Sref, CL_TO, SM, AR_lowspeed, TR);

%--------------------------------------------------------------------------
% performance:
[T_TO_single_engine, ~] = Propulsion(M_TO, alt_TO);
T_TO = T_TO_single_engine*ne; % [N] total takeoff thrust

[S_G, S_TO, BFL] = perf_takeoff(ne, V_stall, V_TO, T_TO, alt_TO, MTOW, Sref, CL_TO, CD_TO);

%--------------------------------------------------------------------------
% S&C:
[CMa, Cl_beta, Cn_beta, Cn_dr, S_VT, l_VT] = stability(M_TO, AR_lowspeed, 24, Sref, b_unswept, TR, CL_TO, SM, 'Takeoff');

%% ========================================================================
% Climb Phase
%--------------------------------------------------------------------------
if M_cr > 1
    alt_climb = alt_sub_cr;              
else
    alt_climb = alt_super_cr;
end

[ROC, gamma_climb, S_climb, dt_climb] = perf_climb(alt_climb, CL_max, W_climb_avg, Sref, b_unswept);

%% ========================================================================
% Cruise Phase
%--------------------------------------------------------------------------
[~, ~, ~, rho_cr, son_cr, ~, ~, ~, ~, ~] = ATMO(alt_cr, 'M');
V_cr = M_cr*son_cr; % [m/s]

CL_cr = W_cruise_start/(Sref*0.5*rho_cr*V_cr^2); % lift coefficient cruise

[R_constH, R_CC, TOF_constH, TOF_CC] = perf_cruise(M_cr, alt_cr, W_cruise_start, W_cruise_end, Sref, SM, AR_swept, TSFC, TR);

%--------------------------------------------------------------------------
% S&C:
[CMa, Cl_beta, Cn_beta, Cn_dr, S_VT, l_VT] = stability(M_cr, AR_swept, 24, Sref, b_swept, TR, CL_cr, SM, 'Cruise');

%% ========================================================================
% Descent Phase
%--------------------------------------------------------------------------
V_descend = 1.5*V_stall; % [m/s] placeholder?
[S_descend, dt_descend] = perf_descent(alt_cr, V_descend);

%% ========================================================================
% Landing Phase
%--------------------------------------------------------------------------
M_land = 0.3; % landing Mach number

[T_land_single_engine, ~] = Propulsion(M_land, alt_land);
T_land = T_land_single_engine*ne; % [N] total landing thrust

[S_land, FAR_land] = perf_land(M_land, alt_land, Sref, AR_lowspeed, W_land, CL_max, T_land, TR, SM);

%--------------------------------------------------------------------------
% S&C:
[CMa, Cl_beta, Cn_beta, Cn_dr, S_VT, l_VT] = stability(M_land, AR_lowspeed, 24, Sref, b_unswept, TR, CL_max, SM, 'Landing');

%% ========================================================================
% Total Performance Summary:
%--------------------------------------------------------------------------
% total range and time of flight:
R_total_1 = S_climb + R_constH + S_descend; % [m]
R_total_2 = S_climb + R_CC + S_descend;     % [m]

dt_total_1 = dt_climb + TOF_constH + dt_descend; % [s]
dt_total_2 = dt_climb + TOF_CC + dt_descend;     % [s]
%--------------------------------------------------------------------------
fprintf('\n\n ============================== Total Range and Time of Flight ============================== \n');
fprintf('\n Total Range (including climb and descent) for Constant Altitude Cruise:');
fprintf('\n -------------------------------------------------------------------------------- ');
fprintf('\n Cruise Altitude: h  = %g [m] = %g [ft] \n', alt_cr, convlength(alt_cr,'m','ft'));
fprintf('\n Range Covered During Climb:        R_climb   = %g [km] = %g [miles]', convlength(S_climb,'m','km'), convlength(S_climb,'m','mi'));
fprintf('\n Constant Altitude Cruise Range:    R_cruise  = %g [km] = %g [miles]', convlength(R_constH,'m','km'), convlength(R_constH,'m','mi'));
fprintf('\n Cruise Climb Range:                R_cruise  = %g [km] = %g [miles]', convlength(R_CC,'m','km'), convlength(R_CC,'m','mi'));
fprintf('\n Range Covered During Descent:      R_descend = %g [km] = %g [miles]', convlength(S_descend,'m','km'), convlength(S_descend,'m','mi'));
fprintf('\n\n Total Range:  R = %g [km] = %g [miles] \n', convlength(R_total_1,'m','km'), convlength(R_total_1,'m','mi'));
fprintf('\n Time to Climb:    dt_climb   = %g [min]', dt_climb/60);
fprintf('\n Time to Cruise:   dt_cruise  = %g [min]', TOF_constH/60);
fprintf('\n Time to Descend:  dt_descend = %g [min]', dt_descend/60);
fprintf('\n\n Total Time of Flight:   dt = %g [min]', dt_total_1/60);
fprintf('\n                         dt = %g [hrs]', dt_total_1/3600);
fprintf('\n\n ===================================================================================== \n');

%% ========================================================================
% Cost analysis:
%--------------------------------------------------------------------------
W_A = 0; % what is this?
[ RTDE_Cost ] = costfunky( weights.W_empty, V_cr, T_max_required, M_max, convforce(MTOW,'N','lbf'), weights.W_fuel, R_total_1*0.000621371, dt_climb/3600, dt_descend/3600, TOF_constH/3600, W_A, ne);
