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
rho = 0.002378;                 % slug/ft3
rho_sl = 0.002378;              % slug/ft3
sigfact = rho./rho_sl;          % Assuming SL takeoff
S_to = 10000;                   % takeoff distance, per Shawn, not me
CL_max_to = 1.8;                % Roskam Part 1, Table 3.1
 
WS_to = transpose(linspace(0, 1000, VectorLength));
TW_to = 37.5.*WS_to./(sigfact.*CL_max_to.*S_to);
%Plot Curve 1
linetype='-';
plot(WS_to, TW_to, ['g' linetype],'LineWidth',3);
axis([WSmin 200 0 TWmax])
hold on


%% Eq 2 Req't: Landing
V_app = 135;                %Approach speed in knots
CL_app = 2.0;               %Need this for AS2

WS_arb = ones(VectorLength,1);
WS_land = WS_arb.*(V_app^2/(17.15^2)*sigfact*CL_app);
TW_land = transpose(linspace(TWmin, 10000, VectorLength));
%Plot Curve 2
linetype='-';
plot(WS_land, TW_land, ['b' linetype],'LineWidth',3);
axis([WSmin 200 0 TWmax])


%% Eq 3 Req't: Second Climb Gradient (FAR 25.121)
CL_max = 1.5;                   %From Roskam Table 3.1
CL_sc = CL_max/1.25^2;
L_D_to = 13.0;                 %Need this for OFW
CGR = 0.027;
N = 4;

WS_sc = transpose(linspace(TWmin, 10000, VectorLength));
TW_arb = ones(VectorLength,1);
TW_sc = TW_arb.*(N/(N-1))*(L_D_to^-1 + CGR);
%Plot Curve 3
linetype='-';
plot(WS_sc, TW_sc, ['c' linetype],'LineWidth',3);
axis([WSmin 200 0 TWmax])


%% Eq 4 Req't: Cruise Matching
e = 0.85;
alt_cr = 18000;
AR_unswept = 10;
M_cr = 0.1:.1:2.0;                                %We are sweeping Mach number

for i = 1: length(M_cr)
[T_cr, a_cr, P_cr, rho_cr] = atmosisa(alt_cr);
rho_cr_eng = rho_cr*0.00194032;                 %Convert density

sweep = acosd(0.7/M_cr(i));                %Subsonic sweep angle as f(M)
AR = AR_unswept.*4.*cosd(sweep).^2;        %Swept AR

V_cr = M_cr(i)*a_cr;                       %Velocity in m/s
V_cr_fts = V_cr*3.28084;                %Velocity in ft/s
q_cr = .5*rho_cr_eng*V_cr_fts^2;       %Dyn pressure for cruise
k = 0.93756;                           %Fraction of cruise weight/gross weight
CD0_sub = 0.022;                       %Need this for OFW
CD0_sup = 0.027;                       %Need this for OFW

%Subsonic Cruise
WS_cr = transpose(linspace(TWmin, 1000, VectorLength));
WS_cr_to = WS_cr./k;
TW_cr_reqd(:,i) = CD0_sub*q_cr./(WS_cr_to) + WS_cr_to/(q_cr*pi*AR*e);

%Supersonic Cruise
WS_cr = transpose(linspace(TWmin, 1000, VectorLength));
WS_cr_to = WS_cr./k;
TW_cr_reqd_sup(:,i) = CD0_sup*q_cr./(WS_cr_to) + WS_cr_to/(q_cr*pi*AR*e);


%% Requirement 5: Max Speed Requirement
K=1/(pi*AR*e);
sigfact2=rho_cr_eng/rho_sl;

%Supersonic Curve
WS_Vmax_sup = transpose(linspace(WSmin, WSmax, VectorLength));       %Sets W/s as "x" in T/W equation
TW_Vmax_sup(:,i) = rho_sl.*V_cr_fts^2.*CD0_sup.*(1./(2.*WS_Vmax_sup))+2*K./(rho_cr_eng.*sigfact.*V_cr_fts^2).*WS_Vmax_sup;    %Eq to plot


%Subsonic Curve
WS_Vmax_sub = transpose(linspace(WSmin, WSmax, VectorLength));       %Sets W/s as "x" in T/W equation
TW_Vmax_sub(:,i) = rho_sl.*V_cr_fts^2.*CD0_sub.*(1./(2.*WS_Vmax_sub))+2*K./(rho_cr_eng.*sigfact.*V_cr_fts^2).*WS_Vmax_sub;    %Eq to plot

end

%Plot Curve for subsonic max speed
linetype='-';
plot(WS_Vmax_sub,TW_Vmax_sub,['b-' linetype],'LineWidth',3)
axis([WSmin WSmax 0 TWmax])

%Plot Curve for supersonic max speed
linetype='-';
plot(WS_Vmax_sup,TW_Vmax_sup,['m-' linetype],'LineWidth',3)
axis([WSmin WSmax 0 TWmax])

%Plot Curve for subsonic cruise
linetype='-';
plot(WS_cr_to, TW_cr_reqd, ['m' linetype],'LineWidth',3);
axis([WSmin 200 0 TWmax])

%Plot Curve for supersonic cruise
linetype='-';
plot(WS_cr_to, TW_cr_reqd_sup, ['b' linetype],'LineWidth',3);
axis([WSmin 200 0 TWmax])

% %Plot Vehicle Markers
% idx = find(abs(WS_land - WS_Vmax_sup) < 1, 1); %// Index of coordinate in array
% py = TW_Vmax_sup(idx);
% px = WS_land(idx);
% designpt = [px,py];
% plot(px,py,'-o','MarkerSize',12, 'MarkerEdgeColor','red','LineWidth',4)
% set(gca,'FontSize',16);
% xlabel('Wing Loading (lb/ft^{2})','FontSize',18);
% ylabel('Thrust Loading','FontSize',18);
% title('Parametric Sizing Chart for SSBJ','FontSize',18);
