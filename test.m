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
%%
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

%Vn_Diagram(MTOW, S, altitude, M_cruise, M_max, CL_max);

[weights, weight_fractions] = Weight_Buildup(passengers, crew, V_cruise, M_cruise, range, SFC, LD_cruise);

%[max_load, min_load] = Wing_Loading(b, convforce(MTOW, 'lbf', 'N'), lambda);
