%{ 
Test script! If you see this, you have git set up correctly
This is an example of how to call a couple functions for the V-n diagram
and weight buildup. Variables listed at the top, function calls below.
%}
clear; clc; close all;

MTOW = 121000; % max takeoff weight, lbf
S = 1350; % wing area, ft^2
altitude = [0, 30000, 51000]; % ft
M_cruise = 1.4;
M_max = 1.5;
CL_max = 2.6;
passengers = 12; % number of passengers
crew = 4; % number of crew members
V_cruise = 924; % MPH
range = 4750; % miles
SFC = 0.9;    % Emperical Placeholder for SFC based on Sadraey Table 4.6 for turbojet
LD_cruise = 7; % Cruise lift/drag from Fig 5.3 Nicolai

Vn_Diagram(MTOW, S, altitude, M_cruise, M_max, CL_max);

[weights, weight_fractions] = Weight_Buildup(MTOW, passengers, crew, V_cruise, M_cruise, range, SFC, LD_cruise);
