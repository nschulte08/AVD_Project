%{
    This script is designed to build the solution space for an OFW SSBJ

% We need LD_to and CD0_sub and CD0_super

===========================================================================
%}
close all; clear; clc;
%% ------------------------------------------------------------------------
% Initialize linear requirements for plotting
VectorLength = 500; % generic length          (for plots)
WSmin = 0.01;        % minimum wing loading   (for plots)
WSmax = 150;         % max wing loading       (for plots)
TWmin = 0.01;        % minimum thrust loading (for plots)
TWmax = 0.5;         % max thrust loading     (for plots)
%--------------------------------------------------------------------------
% other inputs:
e = 0.85;                        % oswald
alt_cr_sub = 13000;              % [m]  subsonic cruise altitude
alt_cr_sub = alt_cr_sub*3.28084; % [ft] subsonic cruise altitude
[~, ~, ~, rho_cr_sub, a_cr_sub, ~, ~, ~, ~, sigfact_cr_sub] = ATMO(alt_cr_sub, 'E'); % 'E' = english units
M_cr_sub = 0.95;
% rho_cr = slug/ft^3
% a_cr = ft/s
% sigfact2 = density ratio = rho/rho_SL
alt_cr_super = 15500;                % [m]  supersonic cruise altitude
alt_cr_super = alt_cr_super*3.28084; % [ft] supersonic cruise altitude
[~, ~, ~, rho_cr_super, a_cr_super, ~, ~, ~, ~, sigfact_cr_super] = ATMO(alt_cr_super, 'E'); % 'E' = english units
M_cr_super = 1.4;

%% ========================================================================
AR_unswept = 8:14; % unswept aspect ratio (sweeping values)

for i = 1:length(AR_unswept)
%% Eq 1 Req't: Takeoff
[~, ~, ~, rho_sl, ~, ~, ~, ~, ~, sigfact_SL] = ATMO(0, 'E'); % 'E' = english units
S_to = 10000;    % [ft] takeoff distance, per Shawn, not me
CL_max_to = 1.8; % max Take off CL , Roskam Part 1, Table 3.1
 
WS_to(:,i) = transpose(linspace(0, 1000, VectorLength)); % generic column vector
TW_to(:,i) = 37.5*WS_to(:,i)/(sigfact_SL*CL_max_to*S_to); % take off thrust loading

%--------------------------------------------------------------------------
%% Eq 2 Req't: Landing
V_app = 135;  % Approach speed in knots
CL_app = 2.0; % Per Roskam Table 3.1

WS_arb = ones(VectorLength,1);
WS_land(:,i) = WS_arb.*(V_app^2/(17.15^2)*sigfact_SL*CL_app);
TW_land(:,i) = transpose(linspace(TWmin, 1000, VectorLength));  % landing thrust loading

%--------------------------------------------------------------------------
%% Eq 3 Req't: Second Climb Gradient (FAR 25.121)
CL_max = 1.5;          %From Roskam Table 3.1
CL_sc = CL_max/1.25^2; % Calculated second climb lift coefficient
L_D_to = 13.0;         % Need this for OFW
CGR = 0.027;           % Climb Gradient = deg/100
N = 4;                 % number of engines

WS_sc(:,i) = transpose(linspace(TWmin, 1000, VectorLength));
TW_arb(:,i) = ones(VectorLength,1);
TW_sc(:,i) = TW_arb(:,i).*(N/(N-1))*(L_D_to^-1 + CGR);  % 2nd climb gradient thrust loading

%--------------------------------------------------------------------------
%% Eq 4 Req't: Cruise Matching
k = 0.93756; % Fraction of cruise weight/gross weight

WS_cr(:,i) = transpose(linspace(TWmin, 1000, VectorLength)); % generic vector (Sets W/S as "x" in T/W equation)
WS_cr_to(:,i) = WS_cr(:,i)./k;                                    % generic vector (Sets W/S as "x" in T/W equation)

%--------------------------------------------------------------------------
%Subsonic Cruise
sweep_sub = acosd(0.7/M_cr_sub);            % Subsonic sweep angle as f(M)
AR_sub = AR_unswept(i)*4*cosd(sweep_sub)^2; % Swept AR
V_cr_sub = M_cr_sub*a_cr_sub;        % [ft/s] Velocity
q_cr_sub = .5*rho_cr_sub*V_cr_sub^2; % [lbf/ft^2] Dyn pressure for cruise

CD0_sub = 0.008;                    % Subsonic CD0 from Concorde analysis

TW_cr_reqd(:,i) = CD0_sub*q_cr_sub./(WS_cr_to(:,i)) + WS_cr_to(:,i)/(q_cr_sub*pi*AR_sub*e);

%--------------------------------------------------------------------------
%Supersonic Cruise
sweep_super = acosd(0.7/M_cr_super);            % Subsonic sweep angle as f(M)
AR_super = AR_unswept(i)*4*cosd(sweep_super)^2; % Swept AR
V_cr_super = M_cr_super*a_cr_super;        % [ft/s] Velocity
q_cr_super = .5*rho_cr_super*V_cr_super^2; % [lbf/ft^2] Dyn pressure for cruise

CD0_super = 0.0104;                 % Supersonic CD0 from Concorde analysis

TW_cr_reqd_super(:,i) = CD0_super*q_cr_super./(WS_cr_to(:,i)) + WS_cr_to(:,i)/(q_cr_super*pi*AR_super*e);

%--------------------------------------------------------------------------
%% Requirement 5: Max Speed Requirement
WS_Vmax = transpose(linspace(WSmin, WSmax, VectorLength)); % generic vector (Sets W/S as "x" in T/W equation)
%--------------------------------------------------------------------------
%Subsonic Curve
K_sub = 1/(pi*AR_sub*e);
TW_Vmax_sub(:,i) = rho_cr_sub.*V_cr_sub^2.*CD0_sub.*(1./(2.*WS_Vmax))+2*K_sub./(rho_cr_sub.*sigfact_cr_sub.*V_cr_sub^2).*WS_Vmax; %Eq to plot

%--------------------------------------------------------------------------
%Supersonic Curve
K_super = 1/(pi*AR_super*e);
TW_Vmax_super(:,i) = rho_cr_super.*V_cr_super^2.*CD0_super.*(1./(2.*WS_Vmax))+2*K_super./(rho_cr_super.*sigfact_cr_super.*V_cr_super^2).*WS_Vmax; % Eq to plot

end
%==========================================================================
%% plot results:
for n = 1:length(AR_unswept)
figure_name = sprintf('Solution Space Plot, AR_unswept = %g',AR_unswept(n));
figure('Name',figure_name,'NumberTitle','off','units','normalized','outerposition',[0 0 1 1]);
hold on
%--------------------------------------------------------------------------
% Take-off
plot(WS_to(:,n), TW_to(:,n), 'g-', 'LineWidth',3);
axis([WSmin, WSmax, 0, TWmax])
%--------------------------------------------------------------------------
% Landing
plot(WS_land(:,n), TW_land(:,n), 'k-','LineWidth',3);
axis([WSmin, WSmax, 0, TWmax])
%--------------------------------------------------------------------------
% Second Climb Gradient (FAR 25.121)
plot(WS_sc(:,n), TW_sc(:,n), 'c-','LineWidth',3);
axis([WSmin, WSmax, 0, TWmax])
%--------------------------------------------------------------------------
% subsonic max speed
plot(WS_Vmax,TW_Vmax_sub(:,n),'b-','LineWidth',3)
axis([WSmin, WSmax, 0, TWmax])
%--------------------------------------------------------------------------
% supersonic max speed
plot(WS_Vmax,TW_Vmax_super(:,n),'m-','LineWidth',3)
axis([WSmin, WSmax, 0, TWmax])
%--------------------------------------------------------------------------
% subsonic cruise
plot(WS_cr_to(:,n), TW_cr_reqd(:,n), 'b--' ,'LineWidth',3);
axis([WSmin, WSmax, 0, TWmax])
%--------------------------------------------------------------------------
% supersonic cruise
plot(WS_cr_to(:,n), TW_cr_reqd_super(:,n), 'm--' ,'LineWidth',3);
axis([WSmin, WSmax, 0, TWmax])
%--------------------------------------------------------------------------
% find design point: (which curves to use was identified manually)
Curve_Fit = polyfit(WS_to(:,n), TW_to(:,n), 1); % get the slope and Y intercept of the straight line (takeoff)
Slope_TO = Curve_Fit(1);
Y_Intercept_TO = Curve_Fit(2);
% plot each "x" of the curved line to get a "y" and subtract it from the real "y":
delta_TW(:,n) = abs(TW_Vmax_super(:,n) - (Slope_TO*WS_Vmax + Y_Intercept_TO)); 
[m,idx] = min(delta_TW(:,n)); % index where difference is minimum
px = WS_Vmax(idx,1);
py = TW_Vmax_super(idx,n);
design_pt(n,:) = [px,py,AR_unswept(n)];
plot(px,py,'o','MarkerSize',12, 'MarkerEdgeColor','red','LineWidth',4)
%--------------------------------------------------------------------------
set(gca,'FontSize',16);
xlabel('Wing Loading (lbf/ft^{2})','FontSize',18);
ylabel('Thrust Loading','FontSize',18);
title_name = sprintf('Parametric Sizing Chart for SSBJ (Unswept AR = %g)',AR_unswept(n));
title(title_name,'FontSize',18);
legend('Take-off','Landing','2nd Climb Gradient','Subsonic Max Speed','Supersonic Max Speed','Subsonic Cruise','Supersonic Cruise','Design Point','Location','best');
grid on
hold off
%--------------------------------------------------------------------------
% save figure to file:
current_folder = pwd;
if ~exist('Figures','dir')
    mkdir('Figures')
end
Figures = sprintf('%s/Figures/', current_folder);

fig_name =  sprintf('%s.fig',figure_name);
jpg_name =  sprintf('%s.jpg',figure_name);

cd(Figures)
    savefig(fig_name);
    saveas(gcf,jpg_name)
    close(gcf)
cd(current_folder);
%--------------------------------------------------------------------------
end
%==========================================================================
%% display design point results:
%design_pt(n,:) = [px,py,AR_unswept(n)];
fprintf('\n\n ========================== Solution Space Results  ========================== \n\n');
for m = 1:length(design_pt(:,1));
fprintf('\n\n ------------------------------------------------------------- \n');
fprintf('\n for an unswept aspect ratio: \n AR_unswept = %g', design_pt(m,3));
fprintf('\n\n Design point: ');
fprintf('\n T/W = %g [lbf/lbf]',   design_pt(m,2));
fprintf('\n W/S = %g  [lbf/ft^2]', design_pt(m,1));
end
fprintf('\n\n ============================================================= \n');
