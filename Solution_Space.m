%{
    This script is designed to build the solution space for an OFW SSBJ

% We need LD_to and CD0_sub and CD0_super
---------------------------------------------------------------------------
Inputs:
AR_unswept:     unspet aspect ratio
CL_max:         max lift coefficient
e:              Oswlad
alt_cr_sub:     subsonic cruise altitude [m]
M_cr_sub:       subsonic cruise Mach number
alt_cr_super:   supersonic cruise altitude [m]
M_cr_super:     supersonic cruise Mach number
ne:             number of engines
---------------------------------------------------------------------------
Outputs:
design_point = [wing_loading, thrust_loading, AR_unswept]

===========================================================================
%}
function [] = Solution_Space(AR_unswept, CL_max, e, alt_cr_sub, M_cr_sub, alt_cr_super, M_cr_super, ne)
%% ------------------------------------------------------------------------
% Initialize linear requirements for plotting
VectorLength = 500; % generic length          (for plots)
WSmin = 0.01;        % minimum wing loading   (for plots)
WSmax = 150;         % max wing loading       (for plots)
TWmin = 0.01;        % minimum thrust loading (for plots)
TWmax = 0.5;         % max thrust loading     (for plots)
%--------------------------------------------------------------------------
% convert to english units:
alt_cr_sub = alt_cr_sub*3.28084; % [ft] subsonic cruise altitude
[~, ~, ~, rho_cr_sub, a_cr_sub, ~, ~, ~, ~, sigfact_cr_sub] = ATMO(alt_cr_sub, 'E'); % 'E' = english units

alt_cr_super = alt_cr_super*3.28084; % [ft] supersonic cruise altitude
[~, ~, ~, rho_cr_super, a_cr_super, ~, ~, ~, ~, sigfact_cr_super] = ATMO(alt_cr_super, 'E'); % 'E' = english units

%% ========================================================================
%% Eq 1 Req't: Takeoff
[~, ~, ~, ~, ~, ~, ~, ~, ~, sigma_SL] = ATMO(0, 'E'); % 'E' = english units
%CL_max = 1.8; % max Take off CL , Roskam Part 1, Table 3.1
WS_to = transpose(linspace(0, 1000, VectorLength)); % generic column vector
%--------------------------------------------------------------------------
S_TO_1 = 4000;    % [ft] takeoff distance
TW_to_1 = 37.5*WS_to/(sigma_SL*CL_max*S_TO_1); % take off thrust loading
%--------------------------------------------------------------------------
S_TO_2 = 5000;    % [ft] takeoff distance
TW_to_2 = 37.5*WS_to/(sigma_SL*CL_max*S_TO_2); % take off thrust loading
%--------------------------------------------------------------------------
S_TO_3 = 6000;    % [ft] takeoff distance
TW_to_3 = 37.5*WS_to/(sigma_SL*CL_max*S_TO_3); % take off thrust loading
%--------------------------------------------------------------------------
S_TO_4 = 7000;    % [ft] takeoff distance
TW_to_4 = 37.5*WS_to/(sigma_SL*CL_max*S_TO_4); % take off thrust loading
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% Eq 2 Req't: Landing
V_app = 135;  % Approach speed in knots
CL_app = 2.0; % Per Roskam Table 3.1
WS_arb = ones(VectorLength,1);
WS_land = WS_arb.*(V_app^2/(17.15^2)*sigma_SL*CL_app);
TW_land = transpose(linspace(TWmin, 1000, VectorLength));  % landing thrust loading

%--------------------------------------------------------------------------
%% Eq 3 Req't: Second Climb Gradient (FAR 25.121)
L_D_to = 13.0;         % Need this for OFW
CGR = 0.0297;           % Climb Gradient = deg/100
WS_sc = transpose(linspace(TWmin, 1000, VectorLength));
TW_arb = ones(VectorLength,1);
TW_sc = TW_arb.*(ne/(ne-1))*(L_D_to^-1 + CGR);  % 2nd climb gradient thrust loading

%--------------------------------------------------------------------------
%% Eq 4 Req't: Cruise Matching
k = 0.5533; % Fraction of cruise weight/gross weight
WS_cr = transpose(linspace(TWmin, 1000, VectorLength));
WS_cr_to = WS_cr./k; % generic vector (Sets W/S as "x" in T/W equation)

%--------------------------------------------------------------------------
%Subsonic Cruise
sweep_sub = acosd(0.7/M_cr_sub);            % Subsonic sweep angle as f(M)
AR_sub = AR_unswept*4*cosd(sweep_sub)^2;    % Swept AR
V_cr_sub = M_cr_sub*a_cr_sub;               % [ft/s] Velocity
q_cr_sub = 0.5*rho_cr_sub*V_cr_sub^2;       % [lbf/ft^2] Dyn pressure for cruise
CD0_sub = 0.008;                            % Subsonic CD0 from Concorde analysis
TW_cr_reqd = CD0_sub*q_cr_sub./(WS_cr_to) + WS_cr_to/(q_cr_sub*pi*AR_sub*e);

%--------------------------------------------------------------------------
%Supersonic Cruise
sweep_super = acosd(0.7/M_cr_super);            % Subsonic sweep angle as f(M)
AR_super = AR_unswept*4*cosd(sweep_super)^2;    % Swept AR
V_cr_super = M_cr_super*a_cr_super;             % [ft/s] Velocity
q_cr_super = 0.5*rho_cr_super*V_cr_super^2;     % [lbf/ft^2] Dyn pressure for cruise
CD0_super = 0.015;                              % Supersonic CD0 from Concorde analysis
TW_cr_reqd_super = CD0_super*q_cr_super./(WS_cr_to) + WS_cr_to/(q_cr_super*pi*AR_super*e);

%--------------------------------------------------------------------------
%% Requirement 5: Max Speed Requirement
WS_Vmax = transpose(linspace(WSmin, WSmax, VectorLength)); % generic vector (Sets W/S as "x" in T/W equation)
%--------------------------------------------------------------------------
%Subsonic Curve
K_sub = 1/(pi*AR_sub*e);
TW_Vmax_sub = rho_cr_sub.*V_cr_sub^2.*CD0_sub.*(1./(2.*WS_Vmax))+2*K_sub./(rho_cr_sub.*sigfact_cr_sub.*V_cr_sub^2).*WS_Vmax; %Eq to plot
%--------------------------------------------------------------------------
%Supersonic Curve
K_super = 1/(pi*AR_super*e);
TW_Vmax_super = rho_cr_super.*V_cr_super^2.*CD0_super.*(1./(2.*WS_Vmax))+2*K_super./(rho_cr_super.*sigfact_cr_super.*V_cr_super^2).*WS_Vmax; % Eq to plot

%% ========================================================================
% plot results:
%--------------------------------------------------------------------------
figure_name = sprintf('Solution Space Plot, AR_unswept = %g',AR_unswept);
figure('Name',figure_name,'NumberTitle','off','units','normalized','outerposition',[0 0 1 1]);
hold on
%--------------------------------------------------------------------------
% Take-off
plot(WS_to, TW_to_1,'Color',[0.0, 0.5, 0.0], 'LineStyle','-','LineWidth',3);
axis([WSmin, WSmax, 0, TWmax])

plot(WS_to, TW_to_2,'Color',[0.0, 0.5, 0.0],'LineStyle','--','LineWidth',3);
axis([WSmin, WSmax, 0, TWmax])

plot(WS_to, TW_to_3,'Color',[0.0, 0.5, 0.0],'LineStyle',':','LineWidth',3);
axis([WSmin, WSmax, 0, TWmax])

plot(WS_to, TW_to_4,'Color',[0.0, 0.5, 0.0],'LineStyle','-.','LineWidth',3);
axis([WSmin, WSmax, 0, TWmax])
%--------------------------------------------------------------------------
% Landing
plot(WS_land, TW_land, 'k-','LineWidth',3);
axis([WSmin, WSmax, 0, TWmax])
%--------------------------------------------------------------------------
% Second Climb Gradient (FAR 25.121)
plot(WS_sc, TW_sc, 'c-','LineWidth',3);
axis([WSmin, WSmax, 0, TWmax])
%--------------------------------------------------------------------------
% subsonic max speed
plot(WS_Vmax,TW_Vmax_sub,'b-','LineWidth',3)
axis([WSmin, WSmax, 0, TWmax])
%--------------------------------------------------------------------------
% supersonic max speed
plot(WS_Vmax,TW_Vmax_super,'m-','LineWidth',3)
axis([WSmin, WSmax, 0, TWmax])
%--------------------------------------------------------------------------
% subsonic cruise
plot(WS_cr_to, TW_cr_reqd, 'b--' ,'LineWidth',3);
axis([WSmin, WSmax, 0, TWmax])
%--------------------------------------------------------------------------
% supersonic cruise
plot(WS_cr_to, TW_cr_reqd_super, 'm--' ,'LineWidth',3);
axis([WSmin, WSmax, 0, TWmax])
%--------------------------------------------------------------------------
set(gca,'FontSize',16);
xlabel('Wing Loading (lbf/ft^{2})','FontSize',18);
ylabel('Thrust Loading','FontSize',18);
title_name = sprintf('Parametric Sizing Chart for SSBJ (Unswept AR = %g)',AR_unswept);
title(title_name,'FontSize',18);
%--------------------------------------------------------------------------
TO_1 = sprintf('S_T_O = %g ft', S_TO_1);
TO_2 = sprintf('S_T_O = %g ft', S_TO_2);
TO_3 = sprintf('S_T_O = %g ft', S_TO_3);
TO_4 = sprintf('S_T_O = %g ft', S_TO_4);

legend(TO_1,TO_2,TO_3,TO_4,'Landing','2nd Climb Gradient','Subsonic Max Speed','Supersonic Max Speed','Subsonic Cruise','Supersonic Cruise','Location','EastOutside');
hold off
%--------------------------------------------------------------------------

end
