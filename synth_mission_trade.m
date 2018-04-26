%{
    Integrated synthesis script for MAE 4351
    Aerospace senior design
%--------------------------------------------------------------------------
Team: ARROW
Team members: 
Shawn McCullough, Ben Holden, Nick Schulte, Rustin Farris, Christian Allen
%--------------------------------------------------------------------------
Last modified: 04/20/2018
% =========================================================================
%}
%clear; clc; close all;
function [OUTPUT] = synth_mission_trade(h, M_cr_sub, M_cr_super)
%% ========================================================================
% Configuration and mission iteration:
config_iter = sprintf('h = %g, M_super = %g', h, M_cr_super); % for saving figures

alt_cr_sub   = 16000; % Subsonic cruise altitude (m)
alt_cr_super = h;     % Supersonic cruise altitude (m)

% AR_unswept = AR_in;                       % Unswept aspect ratio
% TR = TR_in;                               % Wing taper ratio

%% ========================================================================
% Mission inputs
%--------------------------------------------------------------------------
%range_sub    = 9800e3; % Subsonic range, (m)
range_super  = 8800e3; % Supersonic range, (m)

% alt_cr_sub   = 17000;  % Subsonic cruise altitude (m)
% alt_cr_super = 17000;  % Supersonic cruise altitude (m)
% M_cr_super   = 1.6;    % Supersonic cruise Mach number
% M_cr_sub     = 0.7;   % Subsonic cruise Mach number

M_max = 2.8;           % maximum Mach number
num_pass = 12;         % Number of passengers
num_crew = 4;          % Number of crew members
%--------------------------------------------------------------------------
alt_TO = 0;   % takeoff Airport altitude [m]
alt_land = 0; % landing Airport altitude [m]
%--------------------------------------------------------------------------
%[~, TSFC_sub] = Propulsion(M_cr_sub,alt_cr_sub); % TSFC = [1/hr]
[~, TSFC_super] = Propulsion(M_cr_super,alt_cr_super); % TSFC = [1/hr]
%--------------------------------------------------------------------------
% design point: (determined from solution space script)
WingLoading = 30;
ThrustLoading = 0.25;

%% ========================================================================
% Empirical inputs
%--------------------------------------------------------------------------
LD_cruise = 9;  % Cruise lift/drag from Fig 5.3 Nicolai
CL_max = 1.8;   % Placeholder, max CL

%% ========================================================================
% Interdisciplinary inputs (Design inputs) 
%--------------------------------------------------------------------------
% Wing geometry:
AR_unswept = 12;                         % Unswept aspect ratio
%AR_lowspeed = AR_unswept;               % Low speed, unswept AR
TR = 0.5;                               % Wing taper ratio
tmax = 2.3;                             % Maximum thickness, based on AS2 cabin dimensions (m)
tcmax = 0.16;                           % T/c max; variable to iterate
c_r = tmax/tcmax;                       % [m] Root chord = max thickness / tcmax ratio
c_t = TR*c_r;                           % [m] Tip chord
b_unswept = (AR_unswept/2)*(c_r + c_t); % [m] wing span
Sref = b_unswept^2/AR_unswept;          % [m^2] wing area
%--------------------------------------------------------------------------
% Aerodynamics and S&C: 
e = 0.85;                 % Oswald
K = 1/(pi*AR_unswept*e);  % Drag K factor
SM = 0.1;                 % static margin
M_perp = 0.7;             % perp Mach #, variable to iterate(?)
%--------------------------------------------------------------------------
% propulsion:
ne = 4; % number of engines
%T_A_sub = T_A_sub*ne;
%T_A_super = T_A_super*ne;

%% ========================================================================
% Solution Space
%--------------------------------------------------------------------------
Solution_Space(AR_unswept, CL_max, e, alt_cr_sub, M_cr_sub, alt_cr_super, M_cr_super, ne);

fprintf('\n\n ========================== Solution Space Results  ========================== \n');
fprintf('\n The chosen design point is: ');
fprintf('\n T/W  = %g [N] [lbf/lbf] ', ThrustLoading);
fprintf('\n W/S  = %g [N] [lbf/ft^2] ', WingLoading);
fprintf('\n\n ============================================================= \n');

MTOW = convforce(WingLoading*Sref*convlength(1,'m','ft')^2, 'lbf', 'N'); % inital value for max takeoff weight, N

%% ========================================================================
% Weights
%--------------------------------------------------------------------------
[WEIGHTS] = weight_converge(MTOW, alt_cr_super, M_cr_super, range_super, TSFC_super, LD_cruise, num_crew, num_pass, ThrustLoading);
%--------------------------------------------------------------------------
MTOW                    = WEIGHTS.MTOW;
W_empty                 = WEIGHTS.W_empty;
W_fuel                  = WEIGHTS.W_fuel;
W_payload_total         = WEIGHTS.W_payload_total;
W_to_end                = WEIGHTS.W_to_end;
%W_climb_end_super       = WEIGHTS.W_climb_end_super;
W_climb_avg_super       = WEIGHTS.W_climb_avg_super;
W_cruise_start_super    = WEIGHTS.W_cruise_start_super;
W_cruise_end_super      = WEIGHTS.W_cruise_end_super;
W_cruise_avg_super      = WEIGHTS.W_cruise_avg_super;
%W_descend_end_super     = WEIGHTS.W_descend_end_super;
W_descend_avg_super     = WEIGHTS.W_descend_avg_super;
%W_climb_end_sub         = WEIGHTS.W_climb_end_sub;
W_climb_avg_sub         = WEIGHTS.W_climb_avg_sub;
W_cruise_start_sub      = WEIGHTS.W_cruise_start_sub;
W_cruise_end_sub        = WEIGHTS.W_cruise_end_sub;
W_cruise_avg_sub        = WEIGHTS.W_cruise_avg_sub;
%W_descend_end_sub       = WEIGHTS.W_descend_end_sub;
W_descend_avg_sub       = WEIGHTS.W_descend_avg_sub;
W_land                  = WEIGHTS.W_land;
W_cruise_avg_avg        = WEIGHTS.W_cruise_avg_avg;
%--------------------------------------------------------------------------
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
fprintf('\n Empty weight:               W_empty         = %g [N] = %g [lbf] ', W_empty, convforce(W_empty,'N','lbf'));
fprintf('\n Fuel weight:                W_fuel          = %g [N] = %g [lbf] ', W_fuel, convforce(W_fuel,'N','lbf'));
fprintf('\n\n =============================================================================== \n');

%% ========================================================================
% Max thrust required:
T_SL_max_required = ThrustLoading*MTOW; % Required thrust for takeoff, N

%% ========================================================================
% New wing parameters:
%--------------------------------------------------------------------------
S_new_ft = convforce(MTOW,'N','lbf')/WingLoading;  % [ft^2]
S_new = S_new_ft*convlength(1,'ft','m')^2;         % [m^2]
b_new = 2*S_new/(c_r + c_t); % [m] unswept
%new_WS = convforce(MTOW,'N','lbf')/(Sref*convlength(1,'m','ft')^2); % [lbf/ft^2]

fprintf('\n\n ======================================== New Wing ====================================== \n');
fprintf('\n Original wing:  S  = %g [m^2] = %g [ft^2] ', Sref, Sref*convlength(1,'m','ft')^2);
fprintf('\n New wing:       S  = %g [m^2] = %g [ft^2] ', S_new, S_new_ft);
fprintf('\n\n Original wing:  b  = %g [m] = %g [ft]', b_unswept, convlength(b_unswept,'m','ft'));
fprintf(  '\n New wing:       b  = %g [m] = %g [ft]', b_new, convlength(b_new,'m','ft'));
% fprintf('\n\n Original design point:   W/S = %g [N] [lbf/ft^2] ', WingLoading);
% fprintf('\n New design point:        W/S = %g [N] [lbf/ft^2] ', new_WS);
fprintf('\n -------------------------------------------------------------------------------- ');

Sref = S_new; % [m^2]
b_unswept = b_new; % [m]

%% ========================================================================
% vertical tail sizing:
%--------------------------------------------------------------------------
[~, ~, ~, rho_cr_super, son_cr_super, ~, ~, ~, ~, ~] = ATMO(alt_cr_super, 'M');
V_cr_super = M_cr_super*son_cr_super; % [m/s]
%--------------------------------------------------------------------------
% Calculate sweep for cruise:
if M_perp < M_cr_super % to avoid complex numbers
    sweep_deg_cr_super = acosd(M_perp/M_cr_super);      % This has to be limited to 70 ish degrees!
    if sweep_deg_cr_super > 70
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
% size the VT:
[VT, VT_plot] = VT_size(sweep_deg_cr_super, Sref, b_swept_cr_super);

S_VT = VT.S_VT;
l_VT = VT.l_VT;
C_VT = VT.C_VT;
y_VT = VT.y_VT;
TR_VT = VT.TR;
AR_VT = VT.AR;
b_VT = VT.b;
cr_VT = VT.cr;
ct_VT = VT.ct;
cbar_VT = VT.cbar;
Z_bar = VT.Z_bar;
SweepLE_VT = VT.SweepLE_VT;

% =========================================================================
% Plot vertical tail requirements:
%--------------------------------------------------------------------------
%{
figure_name = sprintf('Required Total Vertical Tail Area - Cruise');
figure('Name',figure_name,'NumberTitle','off','units','normalized','outerposition',[0 0 1 1]);
plot(VT_plot(1,:),VT_plot(2,:), 'k', 'LineWidth',3);
xlabel('l_V_T (m)'  ,'FontSize',18);
ylabel('S_V_T (m^2)','FontSize',18);
%xlim([5,b_swept_cr_sub]);
title_string = sprintf('Required Total Vertical Tail Area vs Distance of Vertical Tail mac to cg Location');
title(title_string,'FontSize',18);
grid on
fig_save('Figures', figure_name)
%}
%--------------------------------------------------------------------------

%% ========================================================================
% Display initial design parameters:
%--------------------------------------------------------------------------
fprintf('\n\n =================================== Design Parameters ================================== \n');
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
fprintf('\n Required takeoff thrust:   T  = %g [N] = %g [lbf] ', T_SL_max_required, convforce(T_SL_max_required,'N','lbf'));
fprintf('\n Max takeoff weight:        MTOW  = %g [N] = %g [lbf] ', MTOW, convforce(MTOW,'N','lbf'));
fprintf('\n\n ===================================================================================== \n');
fprintf('\n Vertical tail size estimate:');
fprintf('\n\n Total required vertical tail area:');
fprintf('\n S_VT = %g [m^2]', S_VT);
fprintf('\n\n Location of vertial tail along left wing span:');
fprintf('\n (measured from the center line):');
fprintf('\n y = %g [m]', y_VT);
fprintf('\n\n ===================================================================================== \n');


%% ========================================================================
% operational envelopes:
%--------------------------------------------------------------------------
%{
cruise = [M_cr_sub, alt_cr_sub, M_cr_super, alt_cr_super];
op_envelope(cruise, W_cruise_avg_avg, Sref, SM, b_unswept, TR, CL_max, ne, M_perp)

fig_save_name = sprintf ('Op_Env_config_%s', config_iter);
fig_save('_Results', fig_save_name)
%--------------------------------------------------------------------------
Thrust_required_and_available(cruise, W_cruise_avg_avg, Sref, SM, b_unswept, TR, ne, M_perp)

fig_save_name = sprintf ('Thrust_config_%s', config_iter);
fig_save('_Results', fig_save_name)
%}
%--------------------------------------------------------------------------

%Range_trade(W_cruise_start_super, W_cruise_end_super, Sref, SM, b_unswept, AR_unswept, TR, M_perp)

%--------------------------------------------------------------------------
%% ========================================================================
% V-n diagram and wing loading
%--------------------------------------------------------------------------
%{
altitudes = [0, convlength(alt_cr_sub,'m','ft'), convlength(alt_cr_super,'m','ft')]; % array of key altitudes for V-n diagram (ft)
Vn_Diagram(convforce(MTOW,'N','lbf'), Sref*convlength(1,'m','ft')^2, altitudes, M_cr_sub, M_max, CL_max);
%fig_save('Figures', 'Vn Diagram')
[max_load, min_load] = Wing_Loading(b_unswept, MTOW, TR); 
%fig_save('Figures', 'Wing Loading')
%}
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
[CD_TO, CD0_TO] = aerodynamic_drag(alt_TO, M_TO, Sref, CL_TO, SM, AR_swept_TO, TR, sweep_deg_TO);

%--------------------------------------------------------------------------
% Performance:
[T_TO_single_engine, ~] = Propulsion(M_TO, alt_TO);
T_TO = T_TO_single_engine*ne; % [N] total takeoff thrust

[S_G, S_TO, BFL] = perf_takeoff(ne, V_stall, V_TO, T_TO, alt_TO, MTOW, Sref, CL_TO, CD_TO);

%--------------------------------------------------------------------------
%S&C:
[CMa_TO, Cl_beta_TO, Cn_beta_TO, CM_de_TO, Cl_da_TO, Cn_dr_TO] = stability(M_TO, AR_swept_TO, sweep_deg_TO, Sref, b_swept_TO, TR, CL_TO, SM, S_VT, C_VT, 'Takeoff');

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
[CD_cr_sub, ~] = aerodynamic_drag(alt_cr_sub, M_cr_sub, Sref, CL_cr_sub, SM, AR_swept_cr_sub, TR, sweep_deg_cr_sub);
%--------------------------------------------------------------------------
% propulsion:
T_R_sub = CD_cr_sub*0.5*rho_cr_sub*V_cr_sub^2*Sref; % [N] thrust required

[TSFC_sub, AB_sub, ne] = throttledown(T_R_sub/ne, alt_cr_sub, M_cr_sub, ne);
%--------------------------------------------------------------------------
% Performance:
[R_constH_sub, R_CC_sub, TOF_constH_sub, TOF_CC_sub] = perf_cruise(M_cr_sub, alt_cr_sub, W_cruise_start_sub, W_cruise_end_sub, Sref, TSFC_sub, CL_cr_sub, CD_cr_sub);
fprintf('\n -------------------------------------------------------------------------------- ');
fprintf('\n Subsonic thrust required:   T  = %g [N] = %g [lbf] ', T_R_sub, convforce(T_R_sub,'N','lbf'));
fprintf('\n\n ===================================================================================== \n');

%--------------------------------------------------------------------------
% S&C:
[CMa_cr_sub, Cl_beta_cr_sub, Cn_beta_cr_sub, CM_de_cr_sub, Cl_da_cr_sub, Cn_dr_cr_sub] = stability(M_cr_sub, AR_swept_cr_sub, sweep_deg_cr_sub, Sref, b_swept_cr_sub, TR, CL_cr_sub, SM, S_VT, C_VT, 'Subsonic Cruise');

%% ========================================================================
% Supersonic Cruise Phase
%--------------------------------------------------------------------------
% note: the sweep and density are determined above for VT sizing
%--------------------------------------------------------------------------
% Aero:
CL_cr_super = W_cruise_avg_super/(Sref*0.5*rho_cr_super*V_cr_super^2); % lift coefficient cruise
[CD_cr_super, ~] = aerodynamic_drag(alt_cr_super, M_cr_super, Sref, CL_cr_super, SM, AR_swept_cr_super, TR, sweep_deg_cr_super);
%--------------------------------------------------------------------------
% propulsion:
T_R_super = CD_cr_super*0.5*rho_cr_super*V_cr_super^2*Sref; % [N] thrust required

[TSFC_super, AB_super, ne] = throttledown(T_R_super/ne, alt_cr_super, M_cr_super, ne);
%--------------------------------------------------------------------------
% Performance:
[R_constH_super, R_CC_super, TOF_constH_super, TOF_CC_super] = perf_cruise(M_cr_super, alt_cr_super, W_cruise_start_super, W_cruise_end_super, Sref, TSFC_super, CL_cr_super, CD_cr_super);

fprintf('\n -------------------------------------------------------------------------------- ');
fprintf('\n Supersonic thrust required:   T  = %g [N] = %g [lbf] ', T_R_super, convforce(T_R_super,'N','lbf'));
fprintf('\n\n ===================================================================================== \n');
%--------------------------------------------------------------------------
% S&C:
[CMa_cr_super, Cl_beta_cr_super, Cn_beta_cr_super, CM_de_cr_super, Cl_da_cr_super, Cn_dr_cr_super] = stability(M_cr_super, AR_swept_cr_super, sweep_deg_cr_super, Sref, b_swept_cr_super, TR, CL_cr_super, SM, S_VT, C_VT, 'Supersonic Cruise');

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
% Performance
[S_land, FAR_land, V_TD] = perf_land(alt_land, Sref, AR_swept_Land, W_land, CL_max, TR, SM, sweep_deg_Land);
%--------------------------------------------------------------------------
% S&C:
[CMa_L, Cl_beta_L, Cn_beta_L, CM_de_L, Cl_da_L, Cn_dr_L] = stability(M_Land, AR_swept_Land, sweep_deg_Land, Sref, b_swept_Land, TR, CL_max, SM, S_VT, C_VT, 'Landing');

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
%{
% Vertical tail size:
Take_off = S_VT_TO;
Subsonic_Cruise   = S_VT_cr_sub;
Supersonic_Cruise   = S_VT_cr_super;
Landing  = S_VT_L;
SC_VT = table(Take_off,Subsonic_Cruise,Supersonic_Cruise,Landing);
fprintf('\n ----------------------------------------------------------------------------------------------------------------------------- ');
fprintf('\n Required vertical tail size [m^2]: \n\n');
disp(SC_VT);
%}
%--------------------------------------------------------------------------
% Location of vertical tail:
%{
Take_off = l_VT_TO;
Subsonic_Cruise   = l_VT_cr_sub;
Supersonic_Cruise   = l_VT_cr_super;
Landing  = l_VT_L;
SC_l_VT = table(Take_off,Subsonic_Cruise,Supersonic_Cruise,Landing);
fprintf('\n Corresponding distance (in x-direction) between\n cg location and 1/4 chord of vertical tail mac [m]: \n\n');
disp(SC_l_VT);
%}
% fprintf('\n ----------------------------------------------------------------------------------------------------------------------------- ');
% fprintf('\n ----------------------------------------------------------------------------------------------------------------------------- ');
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
fprintf('\n Cruise Altitude:    h  = %g [m] = %g [ft] \n', alt_cr_sub, convlength(alt_cr_sub,'m','ft'));
fprintf('\n Cruise Mach number: M  = %g \n', M_cr_sub);
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
fprintf('\n Cruise Altitude:    h  = %g [m] = %g [ft] \n', alt_cr_super, convlength(alt_cr_super,'m','ft'));
fprintf('\n Cruise Mach number: M  = %g \n', M_cr_super);
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
% Cost analysis:
%--------------------------------------------------------------------------
PAX = num_pass + num_crew;
[ RTDE_Cost, DOC_Cost ] = costfunky(MTOW, W_empty, W_fuel, V_cr_super, ne, convlength(R_total_super,'m','km'), PAX, T_SL_max_required, M_cr_super);

fprintf('\n\n ================================= Cost Results ================================= \n');
fprintf('\n RDTE cost:    RDTE = %g [$]', RTDE_Cost);
fprintf('\n DOC cost:      DOC = %g [$ per km]', DOC_Cost);
fprintf('\n\n ================================================================================ \n');

%% ========================================================================
% output vector: (for spreadsheet)
%--------------------------------------------------------------------------
if isnan(TSFC_sub) || isnan(TSFC_super)
    OUTPUT = cell(41,1);
else
    OUTPUT = {num_pass; alt_cr_sub; M_cr_sub; alt_cr_super; M_cr_super;...
          AR_unswept; TR; tcmax; tmax; SM; ne;...
          WingLoading; ThrustLoading;...
          Sref; b_unswept; c_r; c_t;...
          MTOW/1000; W_fuel/1000; W_empty/1000; W_land/1000;...
          R_total_sub/1000; dt_total_sub/3600; R_total_super/1000; dt_total_super/3600;...
          V_stall; V_TO; BFL; V_approach; V_TD; FAR_land;...
          sweep_deg_TO; sweep_climb_super(2) ; sweep_deg_cr_sub ;sweep_deg_cr_super ; sweep_descend_super(1); sweep_deg_Land;...
          RTDE_Cost; DOC_Cost;...
          AB_sub; AB_super};
end

%%
close all;
end
