close all; clear all;
%This script is designed to build the solution space for OFW SSBJ
%%======================================================================================%%

% Initialize linear requirements for plotting
VectorLength=1000;
WSmin=0.01;
WSmax=200;
TWmin=0.01;
TWmax=2.0;

%% Eq 1 Req't: Takeoff 
rho = 0.002378;                 %slug/ft3
rho_sl = 0.002378;              %slug/ft3
sigfact = rho./rho_sl;          %Assuming SL takeoff
-S_to = 8000;                    %takeoff distance
-CL_max_to = 1.18;               %Need this for AS2
 
WS_to = transpose(linspace(0, 1000, VectorLength));
TW_to = 37.5.*WS_to./(sigfact.*CL_max_to.*S_to);
%Plot Curve 1
linetype='-';
plot(WS_to, TW_to, ['g' linetype],'LineWidth',3);
axis([WSmin 200 0 TWmax])
hold on


%% Eq 2 Req't: Landing
-V_app = 140;            %Approach speed in knots
-CL_app = 1.2;           %Need this for AS2

WS_arb = ones(VectorLength,1);
WS_land = WS_arb.*(V_app^2/(17.15^2)*sigfact*CL_app);
TW_land = transpose(linspace(TWmin, 10000, VectorLength));
%Plot Curve 2
linetype='-';
plot(WS_land, TW_land, ['b' linetype],'LineWidth',3);
axis([WSmin 200 0 TWmax])
hold on


%% Eq 3 Req't: Second Climb Gradient (FAR 25.121)
-CL_max = 1.18;               %Need this for AS2
CL_sc = CL_max/1.25^2;
-L_D_to = 13.0;               %Need this for AS2
-CGR = 0.027;
-N = 3;

WS_sc = transpose(linspace(TWmin, 10000, VectorLength));
TW_arb = ones(VectorLength,1);
TW_sc = TW_arb.*(N/(N-1))*(L_D_to^-1 + CGR);
%Plot Curve 3
linetype='-';
plot(WS_sc, TW_sc, ['c' linetype],'LineWidth',3);
axis([WSmin 200 0 TWmax])
hold on


%% Eq 4 Req't: Cruise Matching
-AR = 2.76;
-e = 0.8;
-alt_cr = 15000;
[T_cr, a_cr, P_cr, rho_cr] = atmosisa(alt_cr);
rho_cr_eng = rho_cr*0.00194032; %Convert density
-M_cr_sup = 1.4;
-M_cr_sub = 0.95;
V_cr_sub = M_cr_sup*a_cr;               %Velocity in m/s
V_cr_sup = M_cr_sub*a_cr;               %Velocity in m/s
V_cr_fts_sub = V_cr_sub*3.28084;        %Velocity in ft/s
V_cr_fts_sup = V_cr_sup*3.28084;        %Velocity in ft/s
q_sup = .5*rho_cr_eng*V_cr_fts_sup^2;   %Dyn pressure for cruise
q_sub = .5*rho_cr_eng*V_cr_fts_sub^2;   %Dyn pressure for cruise
-k = 0.93756;                    %Fraction of cruise weight/gross weight
-CD0_sub = 0.022;                     %Need this for AS2
-CD0_sup = 0.027;                     %Need this for AS2

%Subsonic Cruise
WS_cr_sub = transpose(linspace(TWmin, 1000, VectorLength));
WS_cr_to = WS_cr_sub./k;
TW_cr_reqd_sub = CD0_sub*q_sub./(WS_cr_to) + WS_cr_to/(q_sub*pi*AR*e);
%Plot Curve 4
linetype='-';
plot(WS_cr_to, TW_cr_reqd_sub, ['m' linetype],'LineWidth',3);
axis([WSmin 200 0 TWmax])
hold on
%Supersonic Cruise
WS_cr = transpose(linspace(TWmin, 1000, VectorLength));
WS_cr_to = WS_cr./k;
TW_cr_reqd_sup = CD0_sup*q_sup./(WS_cr_to) + WS_cr_to/(q_sup*pi*AR*e);
%Plot Curve 4
linetype='-';
plot(WS_cr_to, TW_cr_reqd_sup, ['b' linetype],'LineWidth',3);
axis([WSmin 200 0 TWmax])
hold on

%% Requirement 5: Max Speed Requirement
K=1/(pi*AR*e);
sigfact2=rho_cr_eng/rho_sl;
-CD0_sub = 0.022;            %Need this for AS2
-CD0_sup = 0.027;            %Need this for AS2

%Supersonic Curve
WS_Vmax_sup = transpose(linspace(WSmin, WSmax, VectorLength));       %Sets W/s as "x" in T/W equation
TW_Vmax_sup = rho_sl.*V_cr_fts_sup^2.*CD0_sup.*(1./(2.*WS_Vmax_sup))+2*K./(rho_cr_eng.*sigfact.*V_cr_fts_sup^2).*WS_Vmax_sup;    %Eq to plot
linetype='-';
plot(WS_Vmax_sup,TW_Vmax_sup,['m-' linetype],'LineWidth',3)
axis([WSmin WSmax 0 TWmax])
hold on
% %Subsonic Curve
% WS_Vmax_sub = transpose(linspace(WSmin, WSmax, VectorLength));       %Sets W/s as "x" in T/W equation
% TW_Vmax_sub = rho_sl.*V_cr_fts_sup^2.*CD0_sub.*(1./(2.*WS_Vmax_sub))+2*K./(rho_cr_eng.*sigfact.*V_cr_fts_sub^2).*WS_Vmax_sub;    %Eq to plot
% linetype='-';
% plot(WS_Vmax_sub,TW_Vmax_sub,['b-' linetype],'LineWidth',3)
% axis([WSmin WSmax 0 TWmax])
% hold on


%Plot Vehicle Markers
plot(xas2,yas2,'-*','MarkerSize',10, 'MarkerEdgeColor','red','LineWidth',3)
hold on

idx = find(abs(WS_land - WS_Vmax_sup) < 1, 1); %// Index of coordinate in array
py = TW_Vmax_sup(idx);
px = WS_land(idx);
designpt = [px,py];
plot(px,py,'-o','MarkerSize',12, 'MarkerEdgeColor','red','LineWidth',4)
set(gca,'FontSize',16);
hold on
xlabel('Wing Loading (lb/ft^{2})','FontSize',18);
ylabel('Thrust Loading','FontSize',18);
title('Parametric Sizing Chart for SSBJ','FontSize',18);
