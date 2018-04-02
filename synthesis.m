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
M_max = 2.0;                 % (need to update)
num_pass = 12;               % Number of passengers
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
    %TSFC = 1.0;
end
%% ========================================================================
% Empirical inputs
%--------------------------------------------------------------------------
%TSFC = 0.9;     % [1/hr] Empirical Placeholder for SFC (based on Sadraey Table 4.6 for turbojet = 1.0, turbofan = ???, )
LD_cruise = 9;  % Cruise lift/drag from Fig 5.3 Nicolai
CL_max = 1.8;   % Placeholder, max CL

%% ========================================================================
% Interdisciplinary inputs (Design inputs) 
%--------------------------------------------------------------------------
% wing geometry:
AR_unswept = 8;                         % Unswept aspect ratio
AR_lowspeed = AR_unswept;               % Low speed, unswept AR
TR = 0.3;                               % Wing taper ratio
tmax = 2.3;                             % Maximum thickness, based on AS2 cabin dimensions (m)
tcmax = 0.16;                           % T/c max; variable to iterate
c_r = tmax/tcmax;                       % [m] Root chord = max thickness / tcmax ratio
c_t = TR*c_r;                           % [m] Tip chord
b_unswept = (AR_unswept/2)*(c_r + c_t); % [m] wing span
Sref = b_unswept^2/AR_unswept;          % [m^2] wing area
%--------------------------------------------------------------------------
% aerodynamics and S&C: 
e = 0.85;                   % Oswald
K = 1/(pi*AR_unswept*e);	% Drag K factor
SM = 0.1;                   % static margin
M_perp = 0.7;               % perp Mach #, variable to iterate(?)
%--------------------------------------------------------------------------
% propulsion:
ne = 4; % number of engines

%% ========================================================================
% Solution Space
%--------------------------------------------------------------------------
[design_point] = Solution_Space_OFW_integrated(AR_unswept, CL_max, e, alt_sub_cr, M_cr_sub, alt_super_cr, M_cr_super, ne);

%WingLoading   = design_point(1); % W/S (lbf/ft^2)
%ThrustLoading = design_point(2); % T/W (lbf/lbf)

WingLoading   = 35.18;   % W/S (lbf/ft^2)
ThrustLoading = 0.1718;  % T/W (lbf/lbf)

fprintf('\n -------------------------------------------------------------------------------- ');
fprintf('\n Chosen Design point:');
fprintf('\n T/W  = %g [N] [lbf/lbf] ', ThrustLoading);
fprintf('\n W/S  = %g [N] [lbf/ft^2] ', WingLoading);
fprintf('\n -------------------------------------------------------------------------------- ');

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

W_descend_end = MTOW*wt_frac.WF_to*wt_frac.WF_climb*wt_frac.WF_accel*wt_frac.WF_cruise*wt_frac.WF_des; % Wt at end of descent (N)
W_descend_avg = 0.5*(W_cruise_end + W_descend_end); % Average descent wt (N)

W_land = convforce(weights.W_land,'lbf','N'); % Landing Weight (N)

W_Passengers = convforce(weights.W_payload.Passengers,'lbf','N'); % (N)
W_Luggage = convforce(weights.W_payload.Luggage,'lbf','N');       % (N)
W_Crew = convforce(weights.W_payload.Crew,'lbf','N');             % (N)
W_payload_total = sum([W_Passengers, W_Luggage, W_Crew]); % (N)
W_empty = convforce(weights.W_empty,'lbf','N'); % (N)
W_fuel = convforce(weights.W_fuel,'lbf','N');   % (N)

%% ========================================================================
% Display initial design parameters:
%--------------------------------------------------------------------------
fprintf('\n\n ============================= Initial Design Parameters ============================= \n');
fprintf('\n Required wing area:         S  = %g [m^2] = %g [ft^2] ', Sref, Sref*convlength(1,'m','ft')^2);
fprintf('\n Total unswept wing span:    b  = %g [m] = %g [ft]', b_unswept, convlength(b_unswept,'m','ft'));
fprintf('\n Unswept apect ratio:        AR = %g ', AR_unswept);
fprintf('\n Wign taper ratio:           TR = %g ', TR);
fprintf('\n Maximum wing thickness:     t_max = %g [m] = %g [ft]', tmax, convlength(tmax,'m','ft'));
fprintf('\n Thickness to chord ratio:   t/c = %g ', tcmax);
fprintf('\n Root chord:                 c_r = %g [m] = %g [ft]', c_r, convlength(c_r,'m','ft'));
fprintf('\n Tip chord:                  c_t = %g [m] = %g [ft]', c_t, convlength(c_t,'m','ft'));
fprintf('\n -------------------------------------------------------------------------------- ');
%fprintf('\n Cruise wing sweep:          Lambda = %g [deg] \n', sweep_deg);
%fprintf('\n Effective wing span:        b_eff  = %g [m] = %g [ft]', b_swept, convlength(b_swept,'m','ft'));
%fprintf('\n -------------------------------------------------------------------------------- ');
fprintf('\n Static margin:   SM = %g \n', SM);
fprintf('\n Number of engines:   ne = %g ', ne);
fprintf('\n -------------------------------------------------------------------------------- ');
fprintf('\n Required takeoff thrust:   T  = %g [N] = %g [lbf] ', T_max_required, convforce(T_max_required,'N','lbf'));
fprintf('\n Max takeoff weight:        MTOW  = %g [N] = %g [lbf] ', MTOW, convforce(MTOW,'N','lbf'));
fprintf('\n\n ===================================================================================== \n');

fprintf('\n\n =================================== Weights =================================== \n');
fprintf('\n Max takeoff:                     MTOW  = %g [N] = %g [lbf] ', MTOW,             convforce(MTOW,'N','lbf'));
fprintf('\n Takeoff end:                 W_to_end  = %g [N] = %g [lbf] ', W_to_end,         convforce(W_to_end,'N','lbf'));
fprintf('\n Climb ending:             W_climb_end  = %g [N] = %g [lbf] ', W_climb_end,      convforce(W_climb_end,'N','lbf'));
fprintf('\n Climb average:            W_climb_avg  = %g [N] = %g [lbf] ', W_climb_avg,      convforce(W_climb_avg,'N','lbf'));
fprintf('\n Cruise start:          W_cruise_start  = %g [N] = %g [lbf] ', W_cruise_start,   convforce(W_cruise_start,'N','lbf'));
fprintf('\n Cruise end:              W_cruise_end  = %g [N] = %g [lbf] ', W_cruise_end,     convforce(W_cruise_end,'N','lbf'));
fprintf('\n Cruise average:         W_cruise_avg  = %g [N] = %g [lbf] ', W_cruise_avg,     convforce(W_cruise_avg,'N','lbf'));
fprintf('\n Descent ending:         W_descend_end  = %g [N] = %g [lbf] ', W_descend_end,    convforce(W_descend_end,'N','lbf'));
fprintf('\n Descent average:        W_descend_avg  = %g [N] = %g [lbf] ', W_descend_avg,    convforce(W_descend_avg,'N','lbf'));
fprintf('\n Landing:                       W_land  = %g [N] = %g [lbf] ', W_land,           convforce(W_land,'N','lbf'));
fprintf('\n -------------------------------------------------------------------------------- ');
fprintf('\n Total payload weight:       W_payload       = %g [N] = %g [lbf] ', W_payload_total, convforce(W_payload_total,'N','lbf'));
fprintf('\n Passenger weight:           W_passengers    = %g [N] = %g [lbf] ', W_Passengers, convforce(W_Passengers,'N','lbf'));
fprintf('\n Luggage weight:             W_Luggage       = %g [N] = %g [lbf] ', W_Luggage, convforce(W_Luggage,'N','lbf'));
fprintf('\n Crew weight:                W_Crew          = %g [N] = %g [lbf] ', W_Crew, convforce(W_Crew,'N','lbf'));
fprintf('\n Empty weight:               W_empty         = %g [N] = %g [lbf] ', W_empty, convforce(W_empty,'N','lbf'));
fprintf('\n Fuel weight:                W_fuel          = %g [N] = %g [lbf] ', W_fuel, convforce(W_fuel,'N','lbf'));
fprintf('\n\n =============================================================================== \n');

%% ========================================================================
% operational envelope:
%--------------------------------------------------------------------------
cruise = [M_cr_sub, alt_sub_cr, M_cr_super, alt_super_cr];
op_envelope(cruise, W_cruise_avg, Sref, SM, b_unswept, TR, CL_max, ne)

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
% Calculate sweep for TO:
if M_perp < M_TO % to avoid complex numbers
    sweep_deg_TO = acosd(M_perp/M_TO);      % This has to be limited to 70 ish degrees!
    if sweep_deg_TO > 70                    % Limit the sweep angle
        sweep_deg_TO = 70;
    end
    b_swept_TO = b_unswept*cosd(sweep_deg_TO); % Span at sweep angle [m]
    AR_swept_TO = b_swept_TO^2/Sref;           % Swept aspect ratio
else
    sweep_deg_TO = 0;
	b_swept_TO = b_unswept;   % Span at sweep angle [m]
    AR_swept_TO = AR_unswept; % Swept aspect ratio
end
%--------------------------------------------------------------------------
% aero:
CL_TO = 2*MTOW/(0.5*rho_TO*V_TO^2*Sref);   % N/N
[CD_TO, CD0_TO] = aerofunk_drag_2(alt_TO, M_TO, Sref, CL_TO, SM, AR_swept_TO, TR);

%--------------------------------------------------------------------------
% performance:
[T_TO_single_engine, ~] = Propulsion(M_TO, alt_TO);
T_TO = T_TO_single_engine*ne; % [N] total takeoff thrust

[S_G, S_TO, BFL] = perf_takeoff(ne, V_stall, V_TO, T_TO, alt_TO, MTOW, Sref, CL_TO, CD_TO);

%--------------------------------------------------------------------------
% S&C:
[CMa_TO, Cl_beta_TO, Cn_beta_TO, CM_de_TO, Cl_da_TO, Cn_dr_TO, S_VT_TO, l_VT_TO, VT_plot_TO] = stability(M_TO, AR_swept_TO, sweep_deg_TO, Sref, b_swept_TO, TR, CL_TO, SM, 'Takeoff');

%% ========================================================================
% Climb Phase
%--------------------------------------------------------------------------
if M_cr > 1
    alt_climb = alt_sub_cr;              
else
    alt_climb = alt_super_cr;
end

[ROC, gamma_climb, S_climb, dt_climb, sweep_climb] = perf_climb(alt_climb, CL_max, W_climb_avg, Sref, b_unswept, M_perp);

%% ========================================================================
% Cruise Phase
%--------------------------------------------------------------------------
[~, ~, ~, rho_cr, son_cr, ~, ~, ~, ~, ~] = ATMO(alt_cr, 'M');
V_cr = M_cr*son_cr; % [m/s]
%--------------------------------------------------------------------------
% Calculate sweep for cruise:
if M_perp < M_cr % to avoid complex numbers
    sweep_deg_cr = acosd(M_perp/M_cr);      % This has to be limited to 70 ish degrees!
    if sweep_deg_cr > 70                    % Limit the sweep angle
        sweep_deg_cr = 70;
    end
    b_swept_cr = b_unswept*cosd(sweep_deg_cr); % Span at sweep angle [m]
    AR_swept_cr = b_swept_cr^2/Sref;           % Swept aspect ratio
else
    sweep_deg_cr = 0;
	b_swept_cr = b_unswept;   % Span at sweep angle [m]
    AR_swept_cr = AR_unswept; % Swept aspect ratio
end
%--------------------------------------------------------------------------
% aero:
CL_cr = W_cruise_start/(Sref*0.5*rho_cr*V_cr^2); % lift coefficient cruise
%--------------------------------------------------------------------------
% performance:
[R_constH, R_CC, TOF_constH, TOF_CC] = perf_cruise(M_cr, alt_cr, W_cruise_start, W_cruise_end, Sref, SM, AR_swept_cr, TSFC, TR);
%--------------------------------------------------------------------------
% S&C:
[CMa_cr, Cl_beta_cr, Cn_beta_cr, CM_de_cr, Cl_da_cr, Cn_dr_cr, S_VT_cr, l_VT_cr, VT_plot_cr] = stability(M_cr, AR_swept_cr, sweep_deg_cr, Sref, b_swept_cr, TR, CL_cr, SM, 'Cruise');

%% ========================================================================
% Descent Phase
%--------------------------------------------------------------------------
[S_descend, dt_descend, sweep_descend] = perf_descent(alt_cr, M_perp);

%% ========================================================================
% Landing Phase
%--------------------------------------------------------------------------
[~, ~, ~, rho_land, son_land, ~, ~, ~, ~, ~] = ATMO(0, 'm');
V_stall = sqrt(2*W_land/(rho_land*Sref*CL_max)); % [m/s] stall speed
V_approach = 1.2*V_stall;                        % [m/s] approach velocity
M_Land = V_approach/son_land;                    % landing Mach number
%--------------------------------------------------------------------------
% Calculate sweep for landing:
if M_perp < M_Land % to avoid complex numbers
    sweep_deg_Land = acosd(M_perp/M_Land);      % This has to be limited to 70 ish degrees!
    if sweep_deg_Land > 70                    % Limit the sweep angle
        sweep_deg_Land = 70;
    end
    b_swept_Land = b_unswept*cosd(sweep_deg_Land); % Span at sweep angle [m]
    AR_swept_Land = b_swept_Land^2/Sref;           % Swept aspect ratio
else
    sweep_deg_Land = 0;
	b_swept_Land = b_unswept;   % Span at sweep angle [m]
    AR_swept_Land = AR_unswept; % Swept aspect ratio
end
%--------------------------------------------------------------------------
% propulsion:
[T_land_single_engine, ~] = Propulsion(M_Land, alt_land);
T_land = T_land_single_engine*ne; % [N] total landing thrust

if T_land > 0.1*T_TO % limit landing thrust to 10% of max SL Takeoff thrust
    T_land = 0.1*T_TO; % idle thrust approximately?
end

%--------------------------------------------------------------------------
% performance
[S_land, FAR_land] = perf_land(alt_land, Sref, AR_swept_Land, W_land, CL_max, T_land, TR, SM);
%--------------------------------------------------------------------------
% S&C:
[CMa_L, Cl_beta_L, Cn_beta_L, CM_de_L, Cl_da_L, Cn_dr_L, S_VT_L, l_VT_L, VT_plot_land] = stability(M_Land, AR_swept_Land, sweep_deg_Land, Sref, b_swept_Land, TR, CL_max, SM, 'Landing');

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
% Sweep schedule:
%--------------------------------------------------------------------------
fprintf('\n\n ============================== Sweep Schedule ============================== \n');
fprintf('\n Units: Lambda = [deg], b_eff = [m] \n\n');

ROWNAME = {'Take off';'Bottom of climb';'Top of climb';'Cruise';'Top of descent';'Bottom of descent';'Landing'};

Lambda = [sweep_deg_TO; sweep_climb(1); sweep_climb(2); sweep_deg_cr; sweep_descend(1); sweep_descend(2); sweep_deg_Land];
b_eff  = [b_swept_TO; b_unswept*cosd(sweep_climb(1));   b_unswept*cosd(sweep_climb(2));...
          b_swept_cr; b_unswept*cosd(sweep_descend(1)); b_unswept*cosd(sweep_descend(2)); b_swept_Land];

sweep_schedule = table(Lambda, b_eff,'RowNames',ROWNAME);
disp(sweep_schedule);
fprintf('\n\n ===================================================================================== \n');

%% ========================================================================
% Total S&C Summary:
%--------------------------------------------------------------------------
fprintf('\n\n ================================== S&C Results ================================== \n');
% vertical tail size:
Take_off = S_VT_TO;
Cruise   = S_VT_cr;
Landing  = S_VT_L;
SC_VT = table(Take_off,Cruise,Landing);
fprintf('\n -------------------------------------------------------------------------------- ');
fprintf('\n Required vertical tail size [m^2]: \n\n');
disp(SC_VT);
%--------------------------------------------------------------------------
% location of vertical tail:
Take_off = l_VT_TO;
Cruise   = l_VT_cr;
Landing  = l_VT_L;
SC_l_VT = table(Take_off,Cruise,Landing);
fprintf('\n Corresponding distance (in x-direction) between\n cg location and 1/4 chord of vertical tail mac [m]: \n\n');
disp(SC_l_VT);
fprintf('\n -------------------------------------------------------------------------------- ');
fprintf('\n -------------------------------------------------------------------------------- ');
% longitudinal stability:
ROWNAME = {'CM_alpha [1/rad]';'CM_alpha [1/deg]';'Elevator control power, CM_de [1/rad]';'Elevator control power [1/deg]';};
Take_off = [CMa_TO; CMa_TO*pi/180; CM_de_TO; CM_de_TO*pi/180];
Cruise   = [CMa_cr; CMa_cr*pi/180; CM_de_cr; CM_de_cr*pi/180];
Landing  = [CMa_L; CMa_L*pi/180; CM_de_L; CM_de_L*pi/180];
SC_long = table(Take_off,Cruise,Landing,'RowNames',ROWNAME);
fprintf('\n Longitudinal stability: \n\n');
disp(SC_long);
fprintf('\n -------------------------------------------------------------------------------- ');
fprintf('\n -------------------------------------------------------------------------------- ');
% lateral stability:
ROWNAME = {'Cl_beta [1/rad]';'Cl_betaa [1/deg]';'Aileron control power, Cl_da [1/rad]';'Aileron control power [1/deg]';};
Take_off = [Cl_beta_TO; Cl_beta_TO*pi/180; Cl_da_TO; Cl_da_TO*pi/180];
Cruise   = [Cl_beta_cr; Cl_beta_cr*pi/180; Cl_da_cr; Cl_da_cr*pi/180];
Landing  = [Cl_beta_L; Cl_beta_L*pi/180; Cl_da_L; Cl_da_L*pi/180];
SC_lat = table(Take_off,Cruise,Landing,'RowNames',ROWNAME);
fprintf('\n Lateral stability: \n\n');
disp(SC_lat);
fprintf('\n -------------------------------------------------------------------------------- ');
fprintf('\n -------------------------------------------------------------------------------- ');
% directional stability:
ROWNAME = {'Cn_beta [1/rad]';'Cn_beta [1/deg]';'Rudder control power, Cl_da [1/rad]';'Rudder control power [1/deg]';};
Take_off = [Cn_beta_L; Cn_beta_L*pi/180; Cn_dr_L; Cn_dr_L*pi/180];
Cruise   = [Cn_beta_L; Cn_beta_L*pi/180; Cn_dr_L; Cn_dr_L*pi/180];
Landing  = [Cn_beta_L; Cn_beta_L*pi/180; Cn_dr_L; Cn_dr_L*pi/180];
SC_dir = table(Take_off,Cruise,Landing,'RowNames',ROWNAME);
fprintf('\n Directional stability: \n\n');
disp(SC_dir);
fprintf('\n\n ================================================================================= \n\n');

% =========================================================================
% Plot vertical tail requirements:
%--------------------------------------------------------------------------
if supersonic == 0
    figure_name = sprintf('Required Total Vertical Tail Area - Subsonic Cruise');
elseif supersonic == 1
    figure_name = sprintf('Required Total Vertical Tail Area - Supersonic Cruise');
end
figure('Name',figure_name,'NumberTitle','off','units','normalized','outerposition',[0 0 1 1]);
plot(VT_plot_cr(1,:),VT_plot_cr(2,:), 'k', 'LineWidth',3);
xlabel('l_V_T (m)'  ,'FontSize',18);
ylabel('S_V_T (m^2)','FontSize',18);
xlim([5,b_swept_cr]);
title_string = sprintf('Cruise Required Total Vertical Tail Area vs Distance of Vertical Tail mac to cg Location');
title(title_string,'FontSize',18);
grid on
%fig_save('Figures', figure_name)
%--------------------------------------------------------------------------
figure_name = sprintf('Required Total Vertical Tail Area - TO and Land');
figure('Name',figure_name,'NumberTitle','off','units','normalized','outerposition',[0 0 1 1]);
hold on
plot(VT_plot_TO(1,:),VT_plot_TO(2,:), '--k', 'LineWidth',3);
plot(VT_plot_land(1,:),VT_plot_land(2,:), '--*k', 'LineWidth',3);
hold off
xlabel('l_V_T (m)'  ,'FontSize',18);
ylabel('S_V_T (m^2)','FontSize',18);
xlim([20,b_unswept]);
title_string = sprintf('Required Total Vertical Tail Area vs Distance of Vertical Tail mac to cg Location');
title(title_string,'FontSize',18);
legend('Takeoff','Landing','Location','best');
grid on
%fig_save('Figures', figure_name)

%% ========================================================================
% Cost analysis:
%--------------------------------------------------------------------------
%W_A = 0; % what is this?
%[ RTDE_Cost ] = costfunky( weights.W_empty, V_cr, T_max_required, M_max, convforce(MTOW,'N','lbf'), weights.W_fuel, R_total_1*0.000621371, dt_climb/3600, dt_descend/3600, TOF_constH/3600, W_A, ne);
