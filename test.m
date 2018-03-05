%{ 
Test script! If you see this, you have git set up correctly
This is an example of how to call a couple functions for the V-n diagram
and weight buildup. Variables listed at the top, function calls below.

and some woo hoos.....
%}
clear; clc; close all;
%% ========================================================================
% woooooo hooooooo!!!! 
%--------------------------------------------------------------------------
%{
    my favorite section breaking styles and section comment styles
% -------------------------------------------------------------------------
        see how slick they look!
% -------------------------------------------------------------------------
            woo hoo
% -------------------------------------------------------------------------
                woo hoo
% -------------------------------------------------------------------------
                    woo hoo
% -------------------------------------------------------------------------
%}
% ========================================================================
%% OFW Parameters?
MTOW = 676000; % max takeoff weight, lbf
S = 15700; % wing area, ft^2
altitude = [0, 30000, 51000]; % ft
M_cruise = 1.4;
M_max = 1.5;
CL_max = 2.6;
passengers = 12; % number of passengers
crew = 4; % number of crew members
V_cruise = 923.85; % MPH
range = 5468; % miles (supersonic)
SFC = 1.0;    % Empirical Placeholder for SFC based on Sadraey Table 4.6 for turbojet
LD_cruise = 9; % Cruise lift/drag from Fig 5.3 Nicolai
b = 120; % wingspan, m
lambda = 0.44; % wing taper ratio
K_ng = 1.017; % for pylon-mounted nacelle
N_Lt = 11; % nacelle length, ft
N_w = 4.5; % nacelle width, ft
N_z = 3.75; % ultimate load factor
W_ec = 2.311 * 4500^0.901; % weight of engine and contents, lbf
N_en = 4; % number of engines
S_n =  80; % nacelle wetted area, ft^2

%Vn_Diagram(MTOW, S, altitude, M_cruise, M_max, CL_max);

[weights, weight_fractions] = Weight_Buildup(passengers, crew, V_cruise, M_cruise, range, SFC, LD_cruise);

%[max_load, min_load] = Wing_Loading(b, convforce(MTOW, 'lbf', 'N'), lambda);

%% Concorde parameters
[c_weights, c_weight_fractions] = Weight_Buildup(92, 9, 1354, 2, 4500, 1.195, 7.14);

%% Component weights
W_nacelle = Weight_Component( 1.017, N_Lt, N_w, N_z, W_ec, N_en, S_n );

%% OLD SOLUTION SPACE DESIGN POINT CODE
% find design point: (which curves to use was identified manually)
% Curve_Fit = polyfit(WS_to(:,n), TW_to(:,n), 1); % get the slope and Y intercept of the straight line (takeoff)
% Slope_TO = Curve_Fit(1);
% Y_Intercept_TO = Curve_Fit(2);
% % plot each "x" of the curved line to get a "y" and subtract it from the real "y":
% delta_TW(:,n) = abs(TW_Vmax_super(:,n) - (Slope_TO*WS_Vmax + Y_Intercept_TO)); 
% [m,idx] = min(delta_TW(:,n)); % index where difference is minimum
% px = WS_Vmax(idx,1);
% py = TW_Vmax_super(idx,n);
% design_pt(n,:) = [px,py,AR_unswept(n)];
% plot(px,py,'o','MarkerSize',12, 'MarkerEdgeColor','red','LineWidth',4)
%--------------------------------------------------------------------------
%% display design point results:
% %design_pt(n,:) = [px,py,AR_unswept(n)];
% fprintf('\n\n ========================== Solution Space Results  ========================== \n\n');
% for m = 1:length(design_pt(:,1));
% fprintf('\n\n ------------------------------------------------------------- \n');
% fprintf('\n for an unswept aspect ratio: \n AR_unswept = %g', design_pt(m,3));
% fprintf('\n\n Design point: ');
% fprintf('\n T/W = %g [lbf/lbf]',   design_pt(m,2));
% fprintf('\n W/S = %g  [lbf/ft^2]', design_pt(m,1));
% end
% fprintf('\n\n ============================================================= \n');