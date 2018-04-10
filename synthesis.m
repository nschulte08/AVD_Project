%{
    Integrated synthesis script for MAE 4351
    Aerospace senior design
%--------------------------------------------------------------------------
Team: ARROW
Team members: 
Shawn McCullough, Ben Holden, Nick Schulte, Rustin Farris, Christian Allen
%--------------------------------------------------------------------------
Last modified: 04/05/2018
% =========================================================================
%}
clear; clc; close all;
%% ========================================================================
% Configuration and mission iteration:

config_iter = 'test'; % for saving figures (makes things go much faster)

%% ========================================================================
% Mission inputs
%--------------------------------------------------------------------------
M_cr_sub = 0.78;             % Subsonic cruise Mach number
M_cr_super = 2.0;            % Supersonic cruise Mach number
range_sub = 9800e3;          % Subsonic range, (m)
range_super = 8800e3;        % Supersonic range, (m)
alt_cr_sub = 13000;          % Subsonic cruise altitude (m)
alt_cr_super = 16000;        % Supersonic cruise altitude (m)
M_max = 2.8;                 % (need to update)
num_pass = 12;               % Number of passengers
num_crew = 4;                % Number of crew members
%--------------------------------------------------------------------------
alt_TO = 0;   % takeoff Airport altitude [m]
alt_land = 0; % landing Airport altitude [m]
%--------------------------------------------------------------------------
[~,TSFC_sub] = Propulsion(M_cr_sub,alt_cr_sub); % TSFC = [1/hr]
[~,TSFC_super] = Propulsion(M_cr_super,alt_cr_super); % TSFC = [1/hr]

%% ========================================================================
% Empirical inputs
%--------------------------------------------------------------------------
%TSFC = 0.9;     % [1/hr] Empirical Placeholder for SFC (based on Sadraey Table 4.6 for turbojet = 1.0, turbofan = ???, )
LD_cruise = 9;  % Cruise lift/drag from Fig 5.3 Nicolai
CL_max = 1.8;   % Placeholder, max CL

%% ========================================================================
% Interdisciplinary inputs (Design inputs) 
%--------------------------------------------------------------------------
% Wing geometry:
AR_unswept = 8;                        % Unswept aspect ratio
AR_lowspeed = AR_unswept;               % Low speed, unswept AR
TR = 0.4;                               % Wing taper ratio
tmax = 2.3;                             % Maximum thickness, based on AS2 cabin dimensions (m)
tcmax = 0.16;                           % T/c max; variable to iterate
c_r = tmax/tcmax;                       % [m] Root chord = max thickness / tcmax ratio
c_t = TR*c_r;                           % [m] Tip chord
b_unswept = (AR_unswept/2)*(c_r + c_t); % [m] wing span
Sref = b_unswept^2/AR_unswept;          % [m^2] wing area
%--------------------------------------------------------------------------
% Aerodynamics and S&C: 
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
Solution_Space_OFW_integrated(AR_unswept, CL_max, e, alt_cr_sub, M_cr_sub, alt_cr_super, M_cr_super, ne);

fprintf('\n\n ========================== Solution Space Results  ========================== \n');

fprintf('\n Inspect solution space plot and enter design point: \n');

prompt = ' \n W/S = ';
WingLoading = input(prompt);

prompt = ' \n T/W = ';
ThrustLoading = input(prompt);
%--------------------------------------------------------------------------
% Plot design point on solution space:
hold on
plot(WingLoading,ThrustLoading,'*','MarkerSize',18, 'MarkerEdgeColor','red','LineWidth',3)
legend('Take-off','Landing','2nd Climb Gradient','Subsonic Max Speed','Supersonic Max Speed','Subsonic Cruise','Supersonic Cruise','Design Point','Location','best');
hold off
fig_save_name = sprintf ('Sol_Space_config_%s', config_iter);
fig_save('_Results', fig_save_name)
%--------------------------------------------------------------------------
fprintf('\n\n The chosen design point is: ');
fprintf('\n T/W  = %g [N] [lbf/lbf] ', ThrustLoading);
fprintf('\n W/S  = %g [N] [lbf/ft^2] ', WingLoading);
fprintf('\n\n ============================================================= \n');

MTOW = convforce(WingLoading*Sref*convlength(1,'m','ft')^2, 'lbf', 'N'); % Max takeoff weight, N
T_max_required = ThrustLoading*MTOW; % Required thrust for takeoff, N

%% ========================================================================
% Weights
%--------------------------------------------------------------------------
[~, ~, ~, ~, son_climb_super, ~, ~, ~, ~, ~] = ATMO(alt_cr_super, 'M');
V_cr_super = M_cr_super*son_climb_super; % cruise velocity, (m/s) 
%--------------------------------------------------------------------------
[weights_super, wt_frac_super] = Weight_Buildup(convforce(MTOW,'N','lbf'), num_pass, num_crew, convvel(V_cr_super,'m/s','mph'), M_cr_super, convlength(range_super,'m','mi'), TSFC_super, LD_cruise);

W_to_end = MTOW*wt_frac_super.WF_to;  % Wt at end of TO, start of climb (N)

W_climb_end_super = MTOW*wt_frac_super.WF_to*wt_frac_super.WF_climb*wt_frac_super.WF_accel; % Wt at end of climb, start of cruise (N)
W_climb_avg_super = 0.5*(W_to_end + W_climb_end_super);

W_cruise_start_super = W_climb_end_super; % Wt at beginning of cruise (N)
W_cruise_end_super = MTOW*wt_frac_super.WF_to*wt_frac_super.WF_climb*wt_frac_super.WF_accel*wt_frac_super.WF_cruise; % Wt at end of cruise (N)
W_cruise_avg_super = 0.5*(W_climb_end_super + W_cruise_end_super); % Average cruise wt (N)

W_descend_end_super = MTOW*wt_frac_super.WF_to*wt_frac_super.WF_climb*wt_frac_super.WF_accel*wt_frac_super.WF_cruise*wt_frac_super.WF_des; % Wt at end of descent (N)
W_descend_avg_super = 0.5*(W_cruise_end_super + W_descend_end_super); % Average descent wt (N)

W_land = convforce(weights_super.W_land,'lbf','N'); % Landing Weight (N)
%--------------------------------------------------------------------------
W_Passengers = convforce(weights_super.W_payload.Passengers,'lbf','N'); % (N)
W_Luggage = convforce(weights_super.W_payload.Luggage,'lbf','N');       % (N)
W_Crew = convforce(weights_super.W_payload.Crew,'lbf','N');             % (N)
W_payload_total = sum([W_Passengers, W_Luggage, W_Crew]); % (N)
W_empty = convforce(weights_super.W_empty,'lbf','N'); % (N)
W_fuel = convforce(weights_super.W_fuel,'lbf','N');   % (N)
% =========================================================================
[~, ~, ~, ~, son_climb_sub, ~, ~, ~, ~, ~] = ATMO(alt_cr_sub, 'M');
V_cr_sub = M_cr_sub*son_climb_sub; % cruise velocity, (m/s) 
%--------------------------------------------------------------------------
[~, wt_frac_sub] = Weight_Buildup(convforce(MTOW,'N','lbf'), num_pass, num_crew, convvel(V_cr_sub,'m/s','mph'), M_cr_sub, convlength(range_sub,'m','mi'), TSFC_sub, LD_cruise);

W_climb_end_sub = MTOW*wt_frac_sub.WF_to*wt_frac_sub.WF_climb*wt_frac_sub.WF_accel; % Wt at end of climb, start of cruise (N)
W_climb_avg_sub = 0.5*(W_to_end + W_climb_end_sub);

W_cruise_start_sub = W_climb_end_sub;        % Wt at beginning of cruise (N)
W_cruise_end_sub = W_climb_end_sub - W_fuel; % Wt at end of cruise (N)
W_cruise_avg_sub = 0.5*(W_climb_end_sub + W_cruise_end_sub); % Average cruise wt (N)

W_descend_end_sub = W_cruise_end_sub*wt_frac_sub.WF_des; % Wt at end of cruise (N)
W_descend_avg_sub = 0.5*(W_cruise_end_sub + W_descend_end_sub); % Average descent wt (N)
%--------------------------------------------------------------------------
W_cruise_avg_avg = (W_cruise_avg_sub + W_cruise_avg_super)/2; % used for op_envelope and Ta vs Tr plot

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
fprintf('\n Max takeoff:                MTOW            = %g [N] = %g [lbf] ', MTOW,             convforce(MTOW,'N','lbf'));
fprintf('\n Takeoff end:                W_to_end        = %g [N] = %g [lbf] ', W_to_end,         convforce(W_to_end,'N','lbf'));
fprintf('\n Climb average (sub):        W_climb_avg     = %g [N] = %g [lbf] ', W_climb_avg_sub,      convforce(W_climb_avg_sub,'N','lbf'));
fprintf('\n Climb average (super):      W_climb_avg     = %g [N] = %g [lbf] ', W_climb_avg_super,      convforce(W_climb_avg_super,'N','lbf'));
fprintf('\n -------------------------------------------------------------------------------- ');
fprintf('\n Subsonic Cruise start:      W_cruise_start  = %g [N] = %g [lbf] ', W_cruise_start_sub,   convforce(W_cruise_start_sub,'N','lbf'));
fprintf('\n Subsonic Cruise end:        W_cruise_end    = %g [N] = %g [lbf] ', W_cruise_end_sub,     convforce(W_cruise_end_sub,'N','lbf'));
fprintf('\n Subsonic Cruise average:    W_cruise_avg    = %g [N] = %g [lbf] ', W_cruise_avg_sub,     convforce(W_cruise_avg_sub,'N','lbf'));
fprintf('\n -------------------------------------------------------------------------------- ');
fprintf('\n Supersonic Cruise start:    W_cruise_start  = %g [N] = %g [lbf] ', W_cruise_start_super,   convforce(W_cruise_start_super,'N','lbf'));
fprintf('\n Supersonic Cruise end:      W_cruise_end    = %g [N] = %g [lbf] ', W_cruise_end_super,     convforce(W_cruise_end_super,'N','lbf'));
fprintf('\n Supersonic Cruise average:  W_cruise_avg    = %g [N] = %g [lbf] ', W_cruise_avg_super,     convforce(W_cruise_avg_super,'N','lbf'));
fprintf('\n -------------------------------------------------------------------------------- ');
fprintf('\n Descent average (sub):      W_descend_avg   = %g [N] = %g [lbf] ', W_descend_avg_sub,    convforce(W_descend_avg_sub,'N','lbf'));
fprintf('\n Descent average (super):    W_descend_avg   = %g [N] = %g [lbf] ', W_descend_avg_super,    convforce(W_descend_avg_super,'N','lbf'));
fprintf('\n Landing:                    W_land          = %g [N] = %g [lbf] ', W_land,           convforce(W_land,'N','lbf'));
fprintf('\n -------------------------------------------------------------------------------- ');
fprintf('\n Total payload weight:       W_payload       = %g [N] = %g [lbf] ', W_payload_total, convforce(W_payload_total,'N','lbf'));
fprintf('\n Passenger weight:           W_passengers    = %g [N] = %g [lbf] ', W_Passengers, convforce(W_Passengers,'N','lbf'));
fprintf('\n Luggage weight:             W_Luggage       = %g [N] = %g [lbf] ', W_Luggage, convforce(W_Luggage,'N','lbf'));
fprintf('\n Crew weight:                W_Crew          = %g [N] = %g [lbf] ', W_Crew, convforce(W_Crew,'N','lbf'));
fprintf('\n Empty weight:               W_empty         = %g [N] = %g [lbf] ', W_empty, convforce(W_empty,'N','lbf'));
fprintf('\n Fuel weight:                W_fuel          = %g [N] = %g [lbf] ', W_fuel, convforce(W_fuel,'N','lbf'));
fprintf('\n\n =============================================================================== \n');

%% ========================================================================
% operational envelopes:
%--------------------------------------------------------------------------
cruise = [M_cr_sub, alt_cr_sub, M_cr_super, alt_cr_super];
op_envelope(cruise, W_cruise_avg_avg, Sref, SM, b_unswept, TR, CL_max, ne, M_perp)

fig_save_name = sprintf ('Op_Env_config_%s', config_iter);
fig_save('_Results', fig_save_name)
%--------------------------------------------------------------------------
Thrust_required_and_available(cruise, W_cruise_avg_avg, Sref, SM, b_unswept, TR, ne, M_perp)

fig_save_name = sprintf ('Thrust_config_%s', config_iter);
fig_save('_Results', fig_save_name)
%--------------------------------------------------------------------------

%Range_trade(W_cruise_start_super, W_cruise_end_super, Sref, SM, b_unswept, AR_unswept, TR, M_perp)

%--------------------------------------------------------------------------
%% ========================================================================
% V-n diagram and wing loading
%--------------------------------------------------------------------------
altitudes = [0, convlength(alt_cr_sub,'m','ft'), convlength(alt_cr_super,'m','ft')]; % array of key altitudes for V-n diagram (ft)
Vn_Diagram(convforce(MTOW,'N','lbf'), Sref*convlength(1,'m','ft')^2, altitudes, M_cr_sub, M_max, CL_max);
fig_save('Figures', 'Vn Diagram')
[max_load, min_load] = Wing_Loading(b_unswept, MTOW, TR); 
fig_save('Figures', 'Wing Loading')

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
% Aero:
CL_TO = 2*MTOW/(0.5*rho_TO*V_TO^2*Sref);   % N/N
[CD_TO, CD0_TO] = aerofunk_drag_2(alt_TO, M_TO, Sref, CL_TO, SM, AR_swept_TO, TR);

%--------------------------------------------------------------------------
% Performance:
[T_TO_single_engine, ~] = Propulsion(M_TO, alt_TO);
T_TO = T_TO_single_engine*ne; % [N] total takeoff thrust

[S_G, S_TO, BFL] = perf_takeoff(ne, V_stall, V_TO, T_TO, alt_TO, MTOW, Sref, CL_TO, CD_TO);

%--------------------------------------------------------------------------
% S&C:
[CMa_TO, Cl_beta_TO, Cn_beta_TO, CM_de_TO, Cl_da_TO, Cn_dr_TO, S_VT_TO, l_VT_TO, VT_plot_TO] = stability(M_TO, AR_swept_TO, sweep_deg_TO, Sref, b_swept_TO, TR, CL_TO, SM, 'Takeoff');

%% ========================================================================
% Climb Phase
%--------------------------------------------------------------------------
[ROC_sub, gamma_climb_sub, S_climb_sub, dt_climb_sub, sweep_climb_sub] = perf_climb(alt_cr_sub, CL_max, W_climb_avg_sub, Sref, b_unswept, M_perp);

[ROC_super, gamma_climb_super, S_climb_super, dt_climb_super, sweep_climb_super] = perf_climb(alt_cr_super, CL_max, W_climb_avg_super, Sref, b_unswept, M_perp);

%% ========================================================================
% Subsonic Cruise Phase
%--------------------------------------------------------------------------
[~, ~, ~, rho_cr_sub, son_cr_sub, ~, ~, ~, ~, ~] = ATMO(alt_cr_sub, 'M');
V_cr_sub = M_cr_sub*son_cr_sub; % [m/s]
%--------------------------------------------------------------------------
% Calculate sweep for cruise:
if M_perp < M_cr_sub % to avoid complex numbers
    sweep_deg_cr_sub = acosd(M_perp/M_cr_sub);      % This has to be limited to 70 ish degrees!
    if sweep_deg_cr_sub > 70                    % Limit the sweep angle
        sweep_deg_cr_sub = 70;
    end
    b_swept_cr_sub = b_unswept*cosd(sweep_deg_cr_sub); % Span at sweep angle [m]
    AR_swept_cr_sub = b_swept_cr_sub^2/Sref;           % Swept aspect ratio
else
    sweep_deg_cr_sub = 0;
	b_swept_cr_sub = b_unswept;   % Span at sweep angle [m]
    AR_swept_cr_sub = AR_unswept; % Swept aspect ratio
end
%--------------------------------------------------------------------------
% Aero:
CL_cr_sub = W_cruise_avg_super/(Sref*0.5*rho_cr_sub*V_cr_sub^2); % lift coefficient cruise
%--------------------------------------------------------------------------
% Performance:
[R_constH_sub, R_CC_sub, TOF_constH_sub, TOF_CC_sub] = perf_cruise(M_cr_sub, alt_cr_sub, W_cruise_start_sub, W_cruise_end_sub, Sref, SM, AR_swept_cr_sub, TSFC_sub, TR);
%--------------------------------------------------------------------------
% S&C:
[CMa_cr_sub, Cl_beta_cr_sub, Cn_beta_cr_sub, CM_de_cr_sub, Cl_da_cr_sub, Cn_dr_cr_sub, S_VT_cr_sub, l_VT_cr_sub, VT_plot_cr_sub] = stability(M_cr_sub, AR_swept_cr_sub, sweep_deg_cr_sub, Sref, b_swept_cr_sub, TR, CL_cr_sub, SM, 'Subsonic Cruise');

%% ========================================================================
% Supersonic Cruise Phase
%--------------------------------------------------------------------------
[~, ~, ~, rho_cr_super, son_cr_super, ~, ~, ~, ~, ~] = ATMO(alt_cr_super, 'M');
V_cr_super = M_cr_super*son_cr_super; % [m/s]
%--------------------------------------------------------------------------
% Calculate sweep for cruise:
if M_perp < M_cr_super % to avoid complex numbers
    sweep_deg_cr_super = acosd(M_perp/M_cr_super);      % This has to be limited to 70 ish degrees!
    if sweep_deg_cr_super > 70                    % Limit the sweep angle
        sweep_deg_cr_super = 70;
    end
    b_swept_cr_super = b_unswept*cosd(sweep_deg_cr_super); % Span at sweep angle [m]
    AR_swept_cr_super = b_swept_cr_super^2/Sref;           % Swept aspect ratio
else
    sweep_deg_cr_super = 0;
	b_swept_cr_super = b_unswept;   % Span at sweep angle [m]
    AR_swept_cr_super = AR_unswept; % Swept aspect ratio
end
%--------------------------------------------------------------------------
% Aero:
CL_cr_super = W_cruise_avg_super/(Sref*0.5*rho_cr_super*V_cr_super^2); % lift coefficient cruise
[CD_cr_super, ~] = aerofunk_drag_2(alt_cr_super, M_cr_super, Sref, CL_cr_super, SM, AR_swept_cr_super, TR);
%--------------------------------------------------------------------------
% Performance:
[R_constH_super, R_CC_super, TOF_constH_super, TOF_CC_super] = perf_cruise(M_cr_super, alt_cr_super, W_cruise_start_super, W_cruise_end_super, Sref, SM, AR_swept_cr_super, TSFC_super, TR);
%--------------------------------------------------------------------------
% S&C:
[CMa_cr_super, Cl_beta_cr_super, Cn_beta_cr_super, CM_de_cr_super, Cl_da_cr_super, Cn_dr_cr_super, S_VT_cr_super, l_VT_cr_super, VT_plot_cr_super] = stability(M_cr_super, AR_swept_cr_super, sweep_deg_cr_super, Sref, b_swept_cr_super, TR, CL_cr_super, SM, 'Supersonic Cruise');

%% ========================================================================
% Descent Phase
%--------------------------------------------------------------------------
[S_descend_sub, dt_descend_sub, sweep_descend_sub] = perf_descent(alt_cr_sub, M_perp);

[S_descend_super, dt_descend_super, sweep_descend_super] = perf_descent(alt_cr_super, M_perp);

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
%{
% Propulsion:
[T_land_single_engine, ~] = Propulsion(M_Land, alt_land);
T_land = T_land_single_engine*ne; % [N] total landing thrust

if T_land > 0.1*T_TO % limit landing thrust to 10% of max SL Takeoff thrust
    T_land = 0.1*T_TO; % idle thrust approximately?
end
%}
%--------------------------------------------------------------------------
% Performance
[S_land, FAR_land, V_TD] = perf_land(alt_land, Sref, AR_swept_Land, W_land, CL_max, TR, SM);
%--------------------------------------------------------------------------
% S&C:
[CMa_L, Cl_beta_L, Cn_beta_L, CM_de_L, Cl_da_L, Cn_dr_L, S_VT_L, l_VT_L, VT_plot_land] = stability(M_Land, AR_swept_Land, sweep_deg_Land, Sref, b_swept_Land, TR, CL_max, SM, 'Landing');

%% ========================================================================
% Total Performance Summary:
%--------------------------------------------------------------------------
% Total range and time of flight:
R_total_sub = S_climb_sub   + R_constH_sub   + S_descend_sub;   % [m]
R_total_super = S_climb_super + R_constH_super + S_descend_super; % [m]

dt_total_sub = dt_climb_sub   + TOF_constH_sub   + dt_descend_sub;   % [s]
dt_total_super = dt_climb_super + TOF_constH_super + dt_descend_super; % [s]
%--------------------------------------------------------------------------
fprintf('\n\n ============================== Total Range and Time of Flight for Subsonic Mission ============================== \n');
fprintf('\n Total Range (including climb and descent)');
fprintf('\n -------------------------------------------------------------------------------- ');
fprintf('\n Cruise Altitude: h  = %g [m] = %g [ft] \n', alt_cr_sub, convlength(alt_cr_sub,'m','ft'));
fprintf('\n Range Covered During Climb:        R_climb   = %g [km] = %g [miles]', convlength(S_climb_sub,'m','km'), convlength(S_climb_sub,'m','mi'));
fprintf('\n Constant Altitude Cruise Range:    R_cruise  = %g [km] = %g [miles]', convlength(R_constH_sub,'m','km'), convlength(R_constH_sub,'m','mi'));
%fprintf('\n Cruise Climb Range:                R_cruise  = %g [km] = %g [miles]', convlength(R_CC_sub,'m','km'), convlength(R_CC_sub,'m','mi'));
fprintf('\n Range Covered During Descent:      R_descend = %g [km] = %g [miles]', convlength(S_descend_sub,'m','km'), convlength(S_descend_sub,'m','mi'));
fprintf('\n\n Total Range:  R = %g [km] = %g [miles] \n', convlength(R_total_sub,'m','km'), convlength(R_total_sub,'m','mi'));
fprintf('\n -------------------------------------------------------------------------------- ');
fprintf('\n Time to Climb:    dt_climb   = %g [min]', dt_climb_sub/60);
fprintf('\n Time to Cruise:   dt_cruise  = %g [min]', TOF_constH_sub/60);
fprintf('\n Time to Descend:  dt_descend = %g [min]', dt_descend_sub/60);
fprintf('\n\n Total Time of Flight:   dt = %g [min]', dt_total_sub/60);
fprintf('\n                         dt = %g [hrs]', dt_total_sub/3600);
fprintf('\n\n =================================================================================================================== \n');
%--------------------------------------------------------------------------
fprintf('\n\n ============================== Total Range and Time of Flight for Supersonic Mission ============================== \n');
fprintf('\n Total Range (including climb and descent)');
fprintf('\n -------------------------------------------------------------------------------- ');
fprintf('\n Cruise Altitude: h  = %g [m] = %g [ft] \n', alt_cr_super, convlength(alt_cr_super,'m','ft'));
fprintf('\n Range Covered During Climb:        R_climb   = %g [km] = %g [miles]', convlength(S_climb_super,'m','km'), convlength(S_climb_super,'m','mi'));
fprintf('\n Constant Altitude Cruise Range:    R_cruise  = %g [km] = %g [miles]', convlength(R_constH_super,'m','km'), convlength(R_constH_super,'m','mi'));
%fprintf('\n Cruise Climb Range:                R_cruise  = %g [km] = %g [miles]', convlength(R_CC_super,'m','km'), convlength(R_CC_super,'m','mi'));
fprintf('\n Range Covered During Descent:      R_descend = %g [km] = %g [miles]', convlength(S_descend_super,'m','km'), convlength(S_descend_super,'m','mi'));
fprintf('\n\n Total Range:  R = %g [km] = %g [miles] \n', convlength(R_total_super,'m','km'), convlength(R_total_super,'m','mi'));
fprintf('\n -------------------------------------------------------------------------------- ');
fprintf('\n Time to Climb:    dt_climb   = %g [min]', dt_climb_super/60);
fprintf('\n Time to Cruise:   dt_cruise  = %g [min]', TOF_constH_super/60);
fprintf('\n Time to Descend:  dt_descend = %g [min]', dt_descend_super/60);
fprintf('\n\n Total Time of Flight:   dt = %g [min]', dt_total_super/60);
fprintf('\n                         dt = %g [hrs]', dt_total_super/3600);
fprintf('\n\n =================================================================================================================== \n');

%% ========================================================================
% Sweep schedule:
%--------------------------------------------------------------------------
fprintf('\n\n ============================== Sweep Schedule ============================== \n');
ROWNAME = {'Take off';'Bottom of climb';'Top of climb';'Cruise';'Top of descent';'Bottom of descent';'Landing'};
fprintf('\n Units: Lambda = [deg], b_eff = [m] \n');

fprintf('\n Subsonic mission: \n\n');
Lambda = [sweep_deg_TO; sweep_climb_sub(1); sweep_climb_sub(2); sweep_deg_cr_sub; sweep_descend_sub(1); sweep_descend_sub(2); sweep_deg_Land];
b_eff  = [b_swept_TO; b_unswept*cosd(sweep_climb_sub(1));   b_unswept*cosd(sweep_climb_sub(2));...
          b_swept_cr_sub; b_unswept*cosd(sweep_descend_sub(1)); b_unswept*cosd(sweep_descend_sub(2)); b_swept_Land];

sweep_schedule = table(Lambda, b_eff,'RowNames',ROWNAME);
disp(sweep_schedule);

fprintf('\n -------------------------------------------------------------------------------- ');
fprintf('\n Supersonic mission: \n\n');
Lambda = [sweep_deg_TO; sweep_climb_super(1); sweep_climb_super(2); sweep_deg_cr_super; sweep_descend_super(1); sweep_descend_super(2); sweep_deg_Land];
b_eff  = [b_swept_TO; b_unswept*cosd(sweep_climb_super(1));   b_unswept*cosd(sweep_climb_super(2));...
          b_swept_cr_super; b_unswept*cosd(sweep_descend_super(1)); b_unswept*cosd(sweep_descend_super(2)); b_swept_Land];

sweep_schedule = table(Lambda, b_eff,'RowNames',ROWNAME);
disp(sweep_schedule);
fprintf('\n\n ===================================================================================== \n');

%% ========================================================================
% Total S&C Summary:
%--------------------------------------------------------------------------
fprintf('\n\n ============================================================== S&C Results ============================================================== \n');
% Vertical tail size:
Take_off = S_VT_TO;
Subsonic_Cruise   = S_VT_cr_sub;
Supersonic_Cruise   = S_VT_cr_super;
Landing  = S_VT_L;
SC_VT = table(Take_off,Subsonic_Cruise,Supersonic_Cruise,Landing);
fprintf('\n ----------------------------------------------------------------------------------------------------------------------------- ');
fprintf('\n Required vertical tail size [m^2]: \n\n');
disp(SC_VT);
%--------------------------------------------------------------------------
% Location of vertical tail:
Take_off = l_VT_TO;
Subsonic_Cruise   = l_VT_cr_sub;
Supersonic_Cruise   = l_VT_cr_super;
Landing  = l_VT_L;
SC_l_VT = table(Take_off,Subsonic_Cruise,Supersonic_Cruise,Landing);
fprintf('\n Corresponding distance (in x-direction) between\n cg location and 1/4 chord of vertical tail mac [m]: \n\n');
disp(SC_l_VT);
fprintf('\n ----------------------------------------------------------------------------------------------------------------------------- ');
fprintf('\n ----------------------------------------------------------------------------------------------------------------------------- ');
% Longitudinal stability:
ROWNAME = {'CM_alpha [1/rad]';'CM_alpha [1/deg]';'Elevator control power, CM_de [1/rad]';'Elevator control power [1/deg]';};
Take_off = [CMa_TO; CMa_TO*pi/180; CM_de_TO; CM_de_TO*pi/180];
Subsonic_Cruise   = [CMa_cr_sub; CMa_cr_sub*pi/180; CM_de_cr_sub; CM_de_cr_sub*pi/180];
Supersonic_Cruise   = [CMa_cr_super; CMa_cr_super*pi/180; CM_de_cr_super; CM_de_cr_super*pi/180];
Landing  = [CMa_L; CMa_L*pi/180; CM_de_L; CM_de_L*pi/180];
SC_long = table(Take_off,Subsonic_Cruise,Supersonic_Cruise,Landing,'RowNames',ROWNAME);
fprintf('\n Longitudinal stability: \n\n');
disp(SC_long);
fprintf('\n ----------------------------------------------------------------------------------------------------------------------------- ');
fprintf('\n ----------------------------------------------------------------------------------------------------------------------------- ');
% Lateral stability:
ROWNAME = {'Cl_beta [1/rad]';'Cl_betaa [1/deg]';'Aileron control power, Cl_da [1/rad]';'Aileron control power [1/deg]';};
Take_off = [Cl_beta_TO; Cl_beta_TO*pi/180; Cl_da_TO; Cl_da_TO*pi/180];
Subsonic_Cruise   = [Cl_beta_cr_sub; Cl_beta_cr_sub*pi/180; Cl_da_cr_sub; Cl_da_cr_sub*pi/180];
Supersonic_Cruise   = [Cl_beta_cr_super; Cl_beta_cr_super*pi/180; Cl_da_cr_super; Cl_da_cr_super*pi/180];
Landing  = [Cl_beta_L; Cl_beta_L*pi/180; Cl_da_L; Cl_da_L*pi/180];
SC_lat = table(Take_off,Subsonic_Cruise,Supersonic_Cruise,Landing,'RowNames',ROWNAME);
fprintf('\n Lateral stability: \n\n');
disp(SC_lat);
fprintf('\n ----------------------------------------------------------------------------------------------------------------------------- ');
fprintf('\n ----------------------------------------------------------------------------------------------------------------------------- ');
% Directional stability:
ROWNAME = {'Cn_beta [1/rad]';'Cn_beta [1/deg]';'Rudder control power, Cl_da [1/rad]';'Rudder control power [1/deg]';};
Take_off = [Cn_beta_TO; Cn_beta_TO*pi/180; Cn_dr_TO; Cn_dr_TO*pi/180];
Subsonic_Cruise   = [Cn_beta_cr_sub; Cn_beta_cr_sub*pi/180; Cn_dr_cr_sub; Cn_dr_cr_sub*pi/180];
Supersonic_Cruise   = [Cn_beta_cr_super; Cn_beta_cr_super*pi/180; Cn_dr_cr_super; Cn_dr_cr_super*pi/180];
Landing  = [Cn_beta_L; Cn_beta_L*pi/180; Cn_dr_L; Cn_dr_L*pi/180];
SC_dir = table(Take_off,Subsonic_Cruise,Supersonic_Cruise,Landing,'RowNames',ROWNAME);
fprintf('\n Directional stability: \n\n');
disp(SC_dir);
fprintf('\n\n ========================================================================================================================================= \n\n');

% =========================================================================
% Plot vertical tail requirements:
%--------------------------------------------------------------------------
figure_name = sprintf('Required Total Vertical Tail Area - Cruise');
figure('Name',figure_name,'NumberTitle','off','units','normalized','outerposition',[0 0 1 1]);
hold on
plot(VT_plot_cr_sub(1,:),VT_plot_cr_sub(2,:), 'k', 'LineWidth',3);
plot(VT_plot_cr_super(1,:),VT_plot_cr_super(2,:), '--k', 'LineWidth',3);
hold off
xlabel('l_V_T (m)'  ,'FontSize',18);
ylabel('S_V_T (m^2)','FontSize',18);
xlim([5,b_swept_cr_sub]);
title_string = sprintf('Cruise Required Total Vertical Tail Area vs Distance of Vertical Tail mac to cg Location');
title(title_string,'FontSize',18);
legend('Subsonic','Supersonic','Location','best');
grid on
fig_save('Figures', figure_name)
%--------------------------------------------------------------------------
figure_name = sprintf('Required Total Vertical Tail Area - TO and Land');
figure('Name',figure_name,'NumberTitle','off','units','normalized','outerposition',[0 0 1 1]);
hold on
plot(VT_plot_TO(1,:),VT_plot_TO(2,:), 'k', 'LineWidth',3);
plot(VT_plot_land(1,:),VT_plot_land(2,:), '-+k', 'LineWidth',3);
hold off
xlabel('l_V_T (m)'  ,'FontSize',18);
ylabel('S_V_T (m^2)','FontSize',18);
xlim([20,b_unswept]);
title_string = sprintf('Required Total Vertical Tail Area vs Distance of Vertical Tail mac to cg Location');
title(title_string,'FontSize',18);
legend('Takeoff','Landing','Location','best');
grid on
fig_save('Figures', figure_name)

%% ========================================================================
% Cost analysis:
%--------------------------------------------------------------------------
%W_A = 0; % what is this?
%[ RTDE_Cost ] = costfunky( weights.W_empty, V_cr, T_max_required, M_max, convforce(MTOW,'N','lbf'), weights.W_fuel, R_total_1*0.000621371, dt_climb/3600, dt_descend/3600, TOF_constH/3600, W_A, ne);

%% ========================================================================
% output vector: (for spreadsheet)
%--------------------------------------------------------------------------
OUTPUT = [num_pass; alt_cr_sub; M_cr_sub; alt_cr_super; M_cr_super;...
          AR_unswept; TR; tcmax; tmax; SM; ne;...
          WingLoading; ThrustLoading;...
          Sref; b_unswept; c_r; c_t;...
          MTOW/1000; W_fuel/1000; W_empty/1000; W_land/1000;...
          R_total_sub/1000; dt_total_sub/3600; R_total_super/1000; dt_total_super/3600;...
          V_stall; V_TO; BFL; V_approach; V_TD; FAR_land;...
          sweep_deg_TO; sweep_climb_super(2) ; sweep_deg_cr_sub ;sweep_deg_cr_super ; sweep_descend_super(2); sweep_deg_Land];
